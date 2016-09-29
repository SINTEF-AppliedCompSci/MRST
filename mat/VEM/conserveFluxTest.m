clc; clear; close all;

% w_i are one on one of the faces, 0 at the other. 

G = computeVEMGeometry(voronoiCube(100,[1,1,1]));

rot = @(theta, eta, psi) [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)]...
                        *[cos(eta) 0 sin(eta); 0 1 0; -sin(eta) 0 cos(eta)         ]...
                        *[cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1         ];
K = diag(100*milli*darcy*rand(1,3),0);
R = rot(rand(1)*2*pi, rand(1)*2*pi, rand(1)*2*pi);
K = R'*K*R;
rock.perm = repmat(K([1,2,3,5,6,9]), [G.cells.num, 1]);

fluid      = initSingleFluid('mu' , 1*centi*poise, ...
                             'rho', 1000*kilogram/meter^3);
state = initState(G, [], 0);

p = 1e4;

tol = 1e-6;

bf = boundaryFaces(G);
e = abs(G.faces.centroids(bf,1)-1) < tol;
w = abs(G.faces.centroids(bf,1)) < tol;
bc = addBC([], bf(w), 'pressure', p);
bc = addBC(bc, bf(e), 'pressure', 0);
bc = addBC(bc, bf(~e & ~w), 'flux', 0);

tic

S = computeVirtualIP(G, rock, 2);
stateVEM = incompVEM(state, G, S, fluid, 'bc', bc);

toc

T = computeMultiPointTrans(G, rock);
stateMPFA = incompMPFA(state, G, T, fluid, 'bc', bc);

dd = norm(stateVEM.flux-stateMPFA.flux)/norm(stateMPFA.flux);
fprintf('Flux difference before posprocessing: \t %.2f %%\n', dd*100)

f = G.cells.faces(:,1);
ncf = diff(G.cells.facePos);
fSgn = 1 - 2*(G.faces.neighbors(f,1) ~= rldecode((1:G.cells.num)', ncf,1));

[ii, jj] = blockDiagIndex(ones(G.cells.num,1), ncf);
r = -sparse(ii,jj,1)*(stateVEM.flux(f).*fSgn);

dd = norm(r)/norm(stateVEM.flux);
fprintf('Non-conservativism before: \t\t %.2d\n', dd)

% flux = stateMPFA.flux(f).*fSgn;
% 
% [ii, jj] = blockDiagIndex(ones(G.cells.num,1), ncf);
% fluxMPFA = sparse(ii,jj,1)*flux;
% [ii, jj] = blockDiagIndex(ones(G.cells.num,1), ncf);
% ca = sparse(ii,jj,1)*G.faces.areas(f);

B = zeros(G.cells.num, G.cells.num);


% delta = G.faces.normals(f,:)*K;


ncf = diff(G.cells.facePos);

fn = bsxfun(@times,G.faces.normals(f,:), fSgn./G.faces.areas(f));
delta = (fn*K)';
[ii, jj] = blockDiagIndex(3*ones(numel(f),1), ones(numel(f),1));
delta = sparse(ii,jj,delta(:));

fn = fn';
fn = sparse(ii,jj,fn(:));
delta = fn'*delta;
delta = delta*ones(size(delta,1),1);

ii = f;
jj = (1:numel(f))';
P = sparse(ii, jj,1);
omega = P*delta;

% delta = P.*repmat(delta', G.faces.num,1);
% [~, ~, delta] = find(delta);
% jj = jj(

% delta = prod(delta.*P, 2);
% 
% omega = omega./delta;

% 
for i = 1:G.faces.num
    d = delta(f == i);
    omega(i) = omega(i)/(numel(d)*prod(d));
end
% omega = omega./(2*omega.^2);

% omega = ones(G.faces.num,1);


% fn = fn';
% [ii, jj] = blockDiagIndex(3*ones(G.cells.num,1), ncf);
% fn = sparse(ii,jj,fn(:));
% [ii, jj] = blockDiagIndex(3*ones(G.cells.num,1));
% KK = sparse(ii, jj, repmat(K(:), G.cells.num, 1));
% 
% delta = fn'*KK*fn;

% omega = ones(G.faces.num,1);

for i = 1:G.cells.num
    fi = G.cells.faces(G.cells.facePos(i):G.cells.facePos(i+1)-1);
    for j = 1:G.cells.num
        fj = G.cells.faces(G.cells.facePos(j):G.cells.facePos(j+1)-1);
        ff = fi(ismember(fi, fj));
        ffSgn = 1-2*(G.faces.neighbors(ff,1) ~= i);
        B(i,j) = -sum(1./omega(ff).*G.faces.areas(ff).*ffSgn);
%         B(i,j) = -sum(G.faces.areas(ff));
    end
end

c = G.faces.neighbors(f,:);
c(fSgn == -1,:) = c(fSgn == -1, 2:-1:1);
c = c(:,2);
ncf = diff(G.cells.facePos);
ii = rldecode((1:G.cells.num)', ncf,1);
nz = c~=0;
ii = ii(nz);
c = c(nz);
BB = sparse(ii, c, -1./omega(f(nz)).*G.faces.areas(f(nz)).*fSgn(nz));

[ii, jj] = blockDiagIndex(ones(G.cells.num,1), ncf);
ca = sparse(ii,jj,1)*(1./omega(f).*G.faces.areas(f).*fSgn);

BB = BB + sparse(1:G.cells.num, 1:G.cells.num, -ca);
% ca = sparse(rldecode(1:G.cells.num), ncf, 1), G.faces.areas(f).*fSgn;

B = BB;

beta = B\r;
beta = rldecode(beta, ncf,1);
I = sparse(f, 1:numel(f), 1);
beta = I*beta.*G.faces.areas;

flux = stateVEM.flux - 1./omega.*beta;

dd = norm(flux-stateMPFA.flux)/norm(stateMPFA.flux);
fprintf('Flux difference after posprocessing: \t %.2f %%\n', dd*100)

[ii, jj] = blockDiagIndex(ones(G.cells.num,1), ncf);
r = -sparse(ii,jj,1)*(flux(f).*fSgn);

dd = norm(r)/norm(flux);
fprintf('Non-conservativism after: \t\t %.2d\n', dd)

hold off;

plot(stateVEM.flux, 'o');
hold on
plot(flux, 'd');
plot(stateMPFA.flux, 'sq');

legend('VEM before', 'VEM after', 'MPFA')
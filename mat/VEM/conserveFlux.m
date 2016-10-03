function flux = conserveFlux(state, G, rock, fluid,  bc, src)

if isempty(bc)
    bc = addBC([], boundaryFaces(G), 'flux', 0);
end

f = G.cells.faces(:,1);
ncf = diff(G.cells.facePos);
fSgn = 1 - 2*(G.faces.neighbors(f,1) ~= rldecode((1:G.cells.num)', ncf,1));

[ii, jj] = blockDiagIndex(ones(G.cells.num,1), ncf);
P = sparse(ii, jj, 1);
if ~isempty(src)
    rhs = sparse(src.cell, ones(numel(src.cell),1), src.rate, G.cells.num, 1);
else
    rhs = zeros(G.cells.num,1);
end
r = rhs-P*(state.flux(f).*fSgn);

    if norm(r) > 1e-15

%     tm = totmob(state, fluid);
%     K = bsxfun(@times, permTensor(rock, G.griddim)', tm');
    K = permTensor(rock, G.griddim)';
    [ii, jj] = blockDiagIndex(G.griddim*ones(G.cells.num,1));
        K = sparse(ii, jj, K(:));
    
    fn = bsxfun(@times,G.faces.normals(f,:), fSgn./G.faces.areas(f))';
    [ii, jj] = blockDiagIndex(G.griddim*ones(G.cells.num,1), ncf);
    fnMat = sparse(ii,jj, fn(:));
    delta = fnMat'*K*fnMat;
    delta = spdiags(delta, 0);

    ii = f;
    jj = (1:numel(f))';
    omega = sparse(ii, jj,1)*delta;
    
    for i = 1:G.faces.num
        d = delta(f == i);
        omega(i) = omega(i)/(numel(d)*prod(d));
    end
    
%     omega = ones(G.faces.num,1);

    neu = false(G.faces.num, 1);
    neu(bc.face(strcmp(bc.type, 'flux'))) = true;

    c = G.faces.neighbors(f,:);
    c(fSgn == -1,:) = c(fSgn == -1, 2:-1:1);
%     nz = all(c~=0, 2);
    nz = c(:,2) ~= 0;
%     B = sparse(c(nz,1), c(nz,2), -1./omega(f(nz)).*G.faces.areas(f(nz)).*fSgn(nz));
    B = sparse(c(nz,1), c(nz,2), -G.faces.areas(f(nz))./omega(f(nz)));

    I = ones(numel(f),1);
    I(neu(f)) = 0;
    [ii, jj] = blockDiagIndex(ones(G.cells.num,1), ncf);
%     ca = sparse(ii,jj,1)*(1./omega(f).*G.faces.areas(f).*fSgn);
    ca = sparse(ii,jj,I)*(G.faces.areas(f)./omega(f));
    B = B + spdiags(ca, 0, G.cells.num, G.cells.num);
    
    beta = B\r;
    beta = rldecode(beta, ncf,1);
    beta(neu(f)) = 0;
    
    I = sparse(f,1:numel(f),1);
    beta = I*(beta.*fSgn);%.*sum(I,2);
%     beta = I*(beta);%.*sum(I,2);
    flux = state.flux + G.faces.areas.*beta./omega;
    
    r = rhs-P*(flux(f).*fSgn);
    if norm(r) > 1e-15;
        warning('Could not construct conservative flux field');
    end

    else
        warning('Flux already conservative.');
        flux = state.flux;
    end
    
end

function tm = totmob(state, fluid)

   [mu, ~] = fluid.properties(state);
   s       = fluid.saturation(state);
   kr      = fluid.relperm(s, state);
   
   mob    = bsxfun(@rdivide, kr, mu);
   tm = sum(mob, 2);

end
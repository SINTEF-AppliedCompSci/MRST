mrstModule add vemmech mpsa vem
 

clear all

% Two versions of the O-method is implemented:
%
% 1) The standard version (single continuity point). This breaks down on
% Cartesian grids, and becomes arbitrarily bad when close to Cartesian.
% 2) Generalized version (full displacement continuity enforced weakly).
% This is sort of stable (produce nice-looking results), but fails to
% converge for certain triangular grids, as illustrated below.
%
% A third option is to combine the two options based on the number of cells
% in the interaction region (as a heuristic indicator). This apparently
% works for all grids, but the convergence order is lowered.
%

% number of refinement levels
nref = 4;
% spatial refinement
Nd = 2;

% Switch this on for perturbed grids
pert     = 0;
% Test case number (see definitions below)
testCase = 3;
% Constant (between 0 and 1) in MPSA-W which defines position of continuity point
eta      = 0;

switch testCase
    case 1
        % Cartesian grid, generalized method. This should work
        gridType = 1;
        weakCont = 1;
    case 2
        % Cartesian grid, standard method. Gives decent convergence rate,
        % but local systems are singular
        gridType = 1;
        weakCont = 0;
        eta = 0;
    case 3
        % Triangular grid, 90 degree angles, generalized method. Ok
        gridType = 2;
        weakCont = 1;
    case 4
        % Triangular grid, 90 degree angles, standard method. Breaks down
        gridType = 2;
        weakCont = 0;
        eta = 1/3;
    case 5
        % Triangular grid, 90 degree angles, combined method. Works, may give
        % lower convergence order
        gridType = 2;
        weakCont = 1;
        hybridContPt = 1;
        eta = 1/3;
    case 6
        % Equilateral triangles, generalized method. Does not
        % converge
        gridType = 3;
        weakCont = 1;
    case 7
        % Equilateral triangles, standard method. Works
        gridType = 3;
        weakCont = 0;
        eta = 1/3;
    case 8
        % Equilateral triangles, combined method. Works, may give
        % lower convergence order
        gridType = 3;
        weakCont = 1;
        hybridContPt = 1;
        eta = 1/3;
        
end


% Uniform Lame parameters
mu = 1;
lambda = 1;

% Analytical solution
d1 =  @(x,y) x .* (1-x) .* sin(2 * pi * y);
d2 = @(x,y) sin(2 * pi * x) .* sin(2 * pi * y);
dvec =@(coord) [d1(coord(:,1),coord(:,2)),d2(coord(:,1),coord(:,2))];
%divD = diff(d1,x) + diff(d2,y);
d1dx  =@(x,y) (1-2*x).*sin(2*pi*y);
d1dxdx  =@(x,y) -2.*sin(2*pi*y);
d1dxdy  =@(x,y) (1-2*x).*2*pi.*cos(2*pi*y);d1dydx=d1dxdy;
d2dy  =@(x,y) sin(2*pi*x) .*2*pi.*cos(2*pi*y);
d2dydy  =@(x,y) -sin(2*pi*x) .*(2*pi).^2.*sin(2*pi*y);
d2dydx  =@(x,y) cos(2*pi*x) .*(2*pi).^2.*cos(2*pi*y);d2dxdy=d2dydx;
d1dy  =@(x,y) x.*(1-x).*2*pi.*cos(2*pi*y);
d1dydy  =@(x,y) -x.*(1-x).*(2*pi)^2.*sin(2*pi*y);
d2dx  =@(x,y) cos(2*pi*x).*sin(2*pi*y)*2*pi;
d2dxdx  =@(x,y) -sin(2*pi*x).*sin(2*pi*y)*(2*pi)^2;

divD =@(x,y) d1dx(x,y) + d2dy(x,y);
divDdx =@(x,y) d1dxdx(x,y) + d2dydx(x,y);
divDdy =@(x,y) d1dxdy(x,y) + d2dydy(x,y);    
s11 =@(x,y) 2 * mu * d1dx(x,y) + lambda * divD(x,y);
s11dx=@(x,y) 2*mu *d1dxdx(x,y) + lambda * divDdx(x,y);
s12 =@(x,y)  mu * (d1dy(x,y) + d2dx(x,y));
s12dy =@(x,y) mu * (d1dydy(x,y) + d2dxdy(x,y));
s12dx =@(x,y)  mu * (d1dydx(x,y) + d2dxdx(x,y));
s21 =@(x,y)  mu * (d1dy(x,y) + d2dx(x,y));
s21dx =@(x,y) mu*(d1dydx(x,y) + d2dxdx(x,y));
s21dy =@(x,y)  mu * (d1dydy(x,y) + d2dxdy(x,y));
s22 =@(x,y) 2 * mu * d2dy(x,y) + lambda * divD(x,y);
s22dy =@(x,y) 2 * mu * d2dydy(x,y) + lambda * divDdy(x,y);

d1Full = d1;
d2Full = d2;

s11f = s11;
s21f = s21;
s12f = s12;
s22f = s22;

% Apply divergence
%mrhs1 = inline(simplify(diff(s11,x) + diff(s21,y)));
mrhs1 =@(x,y) s11dx(x,y) + s21dy(x,y);
mrhs2 =@(x,y) s12dx(x,y) + s22dy(x,y);
%mrhs2 = inline(simplify(diff(s12,x) + diff(s22,y))); 
    
for iter1 = 1 : nref
    disp(['Refinement ' num2str(iter1)])
    
    Nx = 2^iter1 * ones(1,Nd);
    G = gridForConvTest(Nx,gridType);
    G = computeVEMGeometry(G);
    G = computeGeometryCalc(G);
    
    Nc = G.cells.num;
    Nf = G.faces.num;
    Nn = G.nodes.num;
    Nd = G.griddim;
    
    isBoundary = any(G.faces.neighbors == 0,2);
    bfn = gridFaceNodes(G,find(isBoundary));
    
    if pert == 1
        xno = G.nodes.coords;
        G.nodes.coords = G.nodes.coords + 0.5 * rand(Nn,Nd) /Nx(1);
        G.nodes.coords(bfn,:) = xno(bfn,:);
        G = computeGeometry(G);
    end
    
    xf = G.faces.centroids;
    xc = G.cells.centroids;
    
    tic;
    muc = mu * ones(Nc,1);
    lambdac = lambda * ones(Nc,1);
    
    constit = shear_normal_stress(Nc, Nd, muc, lambdac, zeros(Nc,1));
    bc = addBC([],find(isBoundary),'pressure',0);
    
    md = MPSA_vectorized(G,constit,'weakCont',weakCont,'bc',bc,'eta',eta,'biot',0);%,0,'bc',bc);
    
    rhsMech = [mrhs1(xc(:,1),xc(:,2)).* G.cells.volumes, mrhs2(xc(:,1),xc(:,2)).* G.cells.volumes];
    san1 = sum([s11f(xf(:,1),xf(:,2)), s21f(xf(:,1),xf(:,2))] .* G.faces.normals,2);
    san2 = sum([s12f(xf(:,1),xf(:,2)), s22f(xf(:,1),xf(:,2))] .* G.faces.normals,2);
    dnum = reshape((md.A \ reshape(rhsMech',[],1)),2,[])';
    toc;
    % Analytical solution
    dex = [d1Full(xc(:,1),xc(:,2)) d2Full(xc(:,1),xc(:,2))];
    
    % Errors in L2 and max-norm
    deL2(iter1) = sqrt(sum(sum(bsxfun(@times,G.cells.volumes.^2,(dex - dnum).^2)))) / sqrt(sum(sum(bsxfun(@times,G.cells.volumes.^2,( dex).^2))));
    dem(iter1) = max(max(abs(dex-dnum)));
    
    stress = md.stress * reshape(dnum',[],1);
    
    s_ex = reshape([san1 , san2]',[],1);
    
    sem(iter1) = max(abs(s_ex - stress))/ max(abs(stress));
    fa = reshape(repmat(G.faces.areas,1,Nd)',[],1);
    seL2(iter1) = sqrt(sum(fa.^2 .* (stress - s_ex).^2)) / sqrt(sum(fa.^2 .* s_ex.^2));
    % make VEM solution
   
    [E,nu] = elasticModuloTransform(lambda,mu,'lam_mu','E_nu');
    %{
    if(G.griddim==2)
      [E,nu]=LMu2ENu_2D(lambda,mu);
    else
      [E,nu]=LMu2ENu_3D(lambda,mu);  
    end
    %} 
    Ev = repmat(E, G.cells.num, 1); 
    nuv = repmat(nu, G.cells.num, 1);
    C=Enu2C(Ev,nuv,G);
    %}
    %{
    C=nan(G.cells.num,numel(constit{1}));
    for i=1:G.cells.num
        C(i,:)=constit{i}(:);
    end
    %}
    % set all boundary to no displacement
    faces= find(any(G.faces.neighbors==0,2));
    inodes = mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces+1)-1);
    nodes = unique(G.faces.nodes(inodes));
    el_bc = struct('disp_bc', struct('nodes', nodes, 'uu', zeros(numel(nodes),G.griddim), 'faces', faces,...
        'uu_face', zeros(numel(nodes),G.griddim), 'mask', true(numel(nodes), G.griddim)), ...
                           'force_bc', []);
    load =@(coord) -[mrhs1(coord(:,1),coord(:,2)), mrhs2(coord(:,1),coord(:,2))];
    
    tic;[uVEM, extra] = VEM_linElast(G, C, el_bc, load);toc;
    deVEM(iter1)=sqrt(sum(sum(bsxfun(@times,(uVEM(G.cells.nodes,:)-dvec(G.nodes.coords(G.cells.nodes,:))).^2,G.weights.cell_nodes),2)))./...
                sqrt(sum(sum(bsxfun(@times,dvec(G.nodes.coords(G.cells.nodes,:)).^2,G.weights.cell_nodes),2)));
    % make CC solution with global interface
    tic;[uu, out] = CC_linElast(G, C, el_bc, load);toc;
    
      
end
%%
log2(deL2(1:end-1)./deL2(2:end))
log2(deVEM(1:end-1)./deVEM(2:end))

return

%%
n=3;
figure(),
subplot(3,2,1)
plotCellData(G,dnum(:,1)),colorbar
subplot(n,2,2)
plotCellData(G,dnum(:,2)),colorbar
subplot(n,2,3)
plotNodeData(G,uVEM(:,1)),colorbar
subplot(n,2,4)
plotNodeData(G,uVEM(:,2)),colorbar
subplot(n,2,5)
plotCellData(G,uu(:,1)),colorbar
subplot(n,2,6)
plotCellData(G,uu(:,2)),colorbar
%%
figure(1),clf,
uu_nn=dvec(G.nodes.coords);
uu_cc=dvec(G.cells.centroids);
val=uVEM-uu_nn;
subplot(2,2,1),
plotNodeData(G,val(:,1));colorbar
subplot(2,2,2),
plotNodeData(G,val(:,2));colorbar
val=dnum-uu_cc;
subplot(2,2,3),
plotCellData(G,val(:,1));colorbar
subplot(2,2,4),
plotCellData(G,val(:,2));colorbar
%%
figure(1),
subplot(2,1,1),
plotCellData(G,mrhs1(G.cells.centroids(:,1),G.cells.centroids(:,2)));colorbar
subplot(2,1,2),
plotCellData(G,mrhs2(G.cells.centroids(:,1),G.cells.centroids(:,2)));colorbar

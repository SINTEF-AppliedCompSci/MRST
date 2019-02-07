mrstModule add vemmech coarsegrid msrsb upr vem

%%

n = 10;
G = pebiGrid(1/n, [2,2]);
G.nodes.coords = G.nodes.coords - 1;

plotGrid(G)
axis equal tight

%%

G = computeGeometry(G);
G = createAugmentedGrid(G);
GG = computeCellDimensions(G);
G = computeCellDimensions2(G);

%%

% Cartesian coarse grid
G_cart = cartGrid([100, 100]);
p_cart = partitionUI(G_cart, [10, 10]);
p_cart = sampleFromBox(G, reshape(p_cart, G_cart.cartDims));
GC = generateCoarseGrid(G, p_cart);
GC = coarsenGeometry(GC);
GC = addCoarseCenterPoints(GC);

%%

GC = coarsenCellDimensions(GC);

%%

kMax = 10;

ic = any(GC.faces.neighbors == 2, 2);

cub = CoarseGrid2DCubature(GC, kMax, ic);

%%

tol = 1e-9;
for k = 1:kMax
    
    if k == 1
        stndrdth = 'st';
    elseif k == 2
        stndrdth = 'nd';
    elseif k == 3
        stndrdth = 'rd';
    else
        stndrdth = 'th';
    end
    
    basis  = dgBasis(GC.griddim, k, 'legendre');
    nDof   = basis.nDof;
    sol    = zeros(nDof,1);
    sol(1) = sum(GC.cells.volumes);

    % Test volume cubatures
    [W, x] = cub.getCubature((1:GC.cells.num)', 'volume');    
    I = zeros(nDof,1);
    for dofNo = 1:nDof
        I(dofNo) = sum(W*basis.psi{dofNo}(x));
    end
    assert(all(abs(I - sol) < tol), ...
        [num2str(k), stndrdth, ' order ', class(cub), ' failed']);
    fprintf([num2str(k), stndrdth, ' order ', class(cub), ' successfull\n']);

end
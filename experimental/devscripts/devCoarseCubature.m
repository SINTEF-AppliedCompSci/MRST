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

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>

clear all
close all

N = 20;
G = cartGrid([N, N], [1 , 1]);
angle = 10/180*pi;
% angle = 0;
[G, bcfaces] = rotateGrid(G, angle);
G = computeGeometry(G);

dim = G.griddim;
nc = G.cells.num;

lambda = ones(nc, 1);
mu     = ones(nc, 1);

prop = struct('lambda', lambda, ...
              'mu', mu);

useVirtual = true;
[tbls, mappings] = setupMpsaStandardTables(G, 'useVirtual', useVirtual);

extfaces{1} = bcfaces.ymin;
n = numel(extfaces{1});
linforms{1} = repmat([sin(angle), -cos(angle)], n, 1);
linformvals{1} = zeros(n, 1);

extfaces{2} = bcfaces.xmin;
n = numel(extfaces{2});
linforms{2} = repmat([cos(angle), sin(angle)], n, 1);
linformvals{2} = zeros(n, 1);

nodefacevectbl = tbls.nodefacevectbl;
vectbl         = tbls.vectbl;
facevectbl     = tbls.facevectbl;
facetbl        = tbls.facetbl;

extfacetbl.faces = bcfaces.ymax;
extfacetbl = IndexArray(extfacetbl);

if useVirtual
    face_from_extface = extfacetbl.get('faces');
end

extfacevectbl = crossIndexArray(extfacetbl, vectbl, {}, 'optpureproduct', true, 'virtual', useVirtual);

normals = G.faces.normals;
normals = reshape(normals', [], 1);

map = TensorMap();
map.fromTbl  = facevectbl;
map.toTbl    = extfacevectbl;
map.mergefds = {'faces', 'vec'};

if useVirtual

    map.pivottbl = extfacevectbl;
    [vec, i] = ind2sub([vectbl.num, extfacetbl.num], (1 : extfacevectbl.num)');
    map.dispind1 = sub2ind([vectbl.num, facetbl.num], vec, face_from_extface(i));
    map.dispind2 = (1 : extfacevectbl.num)';
    map.issetup = true;
    
else
    
    map = map.setup();
    
end

extFacetNormals = map.eval(normals);

map = TensorMap();
map.fromTbl  = extfacevectbl;
map.toTbl    = nodefacevectbl;
map.mergefds = {'faces', 'vec'};

if useVirtual
    
    nodefacetbl = tbls.nodefacetbl;
    [extnodefacetbl, indstruct] = crossIndexArray(extfacetbl, nodefacetbl, {'faces'});
    extface_from_extnodeface  = indstruct{1}.inds;
    nodeface_from_extnodeface = indstruct{2}.inds;
    extnodefacevectbl = crossIndexArray(extnodefacetbl, vectbl, {}, 'optpureproduct', true, 'virtual', useVirtual);
    
    map.pivottbl = extnodefacevectbl;
    [vec, i] = ind2sub([vectbl.num, extnodefacetbl.num], (1 : extnodefacevectbl.num)');
    map.dispind1 = sub2ind([vectbl.num, extfacetbl.num], vec, extface_from_extnodeface(i));
    map.dispind2 = sub2ind([vectbl.num, nodefacetbl.num], vec, nodeface_from_extnodeface(i));
    map.issetup = true;
    
else
    map = map.setup();
end

extforce = map.eval(-extFacetNormals);

cellvectbl = tbls.cellvectbl;
force = zeros(cellvectbl.num, 1);

bc.extfaces    = vertcat(extfaces{:});
bc.linform     = vertcat(linforms{:});
bc.linformvals = vertcat(linformvals{:});

bc = setupFaceBC2(bc, G, tbls, mappings, 'useVirtual', useVirtual);

loadstruct.bc = bc;
loadstruct.extforce = extforce;
loadstruct.force = force;

eta = 1/3;
bcetazero = false;

assembly = assembleMPSA2(G, prop, loadstruct, eta, tbls, mappings, 'bcetazero', ...
                         bcetazero, 'extraoutput', true, 'useVirtual', useVirtual);

B   = assembly.B  ;
rhs = assembly.rhs;

sol = B\rhs;

% displacement values at cell centers.
n = cellvectbl.num;

ucell = sol(1 : n);
lagmult = sol(n + 1 : end);

unode = computeNodeDisp2(ucell, tbls, mappings, 'useVirtual', useVirtual);

dim = G.griddim;
ucell = reshape(ucell, dim, [])';
unode = reshape(unode, dim, [])';

figure(1)
plotGrid(G, 'facecolor', 'none', 'edgecolor', 'blue');
plotGridDeformed(G, unode, 'facecolor', 'none', 'edgecolor', 'red');
axis equal

figure(2)
plotCellData(G, ucell(:, 1));
title('disp x direction');
colorbar
axis equal

figure(3)
plotCellData(G, ucell(:, 2));
title('disp y direction');
colorbar
axis equal

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2020 University of Bergen and SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of the MPSA-W module for the MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% The MPSA-W module is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% The MPSA-W module is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with the MPSA-W module.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>


%% Example: three uses of simpleGrdecl
grdecl = simpleGrdecl([20, 20, 5]);
G = processGRDECL(grdecl);
clf, plotGrid(G,'FaceColor',[.8 .8 .8]), view(3), axis tight off

%%
grdecl = simpleGrdecl([20, 20, 5], @(x) 0.05 * (sin(2*pi*x) - 1.5));
G = processGRDECL(grdecl);
clf, plotGrid(G,'FaceColor',[.8 .8 .8]), view(3), axis tight off

%%
grdecl = simpleGrdecl([20, 20, 5], @(x) 0.25*(x-0.5), 'flat', true);
G = processGRDECL(grdecl);
clf, plotGrid(G,'FaceColor',[.8 .8 .8]), view(3), axis tight off

%% Example: makeModel3
grdecl = makeModel3([30,20,5]);
G = processGRDECL(grdecl);
clf, plotGrid(G,'FaceColor',[.8 .8 .8]), view(130,30), axis tight off
plotFaces(G,find(G.faces.tag>0),'FaceColor',[.6 .6 .6]);

%% Example: extrudedTriangleGrid
G = extrudedTriangleGrid(50);
clf, plotGrid(G,'FaceColor',[.8 .8 .8]), view(-30,60), axis tight off

%% 
G = extrudedTriangleGrid(50, true);
clf, plotGrid(G,'FaceColor',[.8 .8 .8]), view(-30,60), axis tight off

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

%% Simple examples
% In this example we will go through a set of elementary examples that
% demonstrate how to generate petrophysical properties in MRST


%% Homogeneous, isotropic properties
G = cartGrid([1000 100]);
rock = makeRock(G, 200*milli*darcy,.2);
plotCellData(G, rock.poro,'EdgeColor','w'); 

%% Homogeneous, anisotropic properties
rock.perm = repmat( [100 100 10].*milli*darcy, [G.cells.num,1]);
clf
plotCellData(G, rock.perm./max(rock.perm(:)), 'edgeColor', 'w');

%% Carman-Kozeny
% First in 2D
G = cartGrid([50 20]);
p = gaussianField(G.cartDims, [0.2 0.4], [11 3], 2.5);
K = p.^3.*(1e-5)^2./(0.81*72*(1-p).^2);
rock = makeRock(G, K(:), p(:));

clf
plotCellData(G,rock.poro);
colorbar('horiz'); axis equal tight;

%%
% Then in 3D
G = cartGrid([50 20 10]);
p = gaussianField(G.cartDims, [0.2 0.4], [5 3 1], 2.5);
K = p.^3.*(1e-5)^2./(0.81*72*(1-p).^2);
rock = makeRock(G, K(:), p(:));

clf
plotCellData(G,convertTo(rock.perm,milli*darcy));
colorbar('horiz'); axis equal tight; view(3);

%% Lognormal, layered in 3D
G = processGRDECL(simpleGrdecl([50 30 10], 0.12));
K = logNormLayers(G.cartDims, [100 400 50 350], ...
   'indices', [1 2 5 7 11]);

clf
plotCellData(G,log10(K),'EdgeColor','k'); view(45,30);
axis tight off, set(gca,'DataAspect',[0.5 1 1])
h=colorbar('horiz'); ticks=25*2.^(0:5);
set(h,'XTick',log10(ticks),'XTickLabel',ticks);

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

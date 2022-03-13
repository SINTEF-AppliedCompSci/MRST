%% Demonstrate lack of convergence for the TPFA scheme
% We consider a homogeneous, rectangular reservoir with a symmetric well
% pattern consisting of one injector and two producers. Because of the
% symmetry, the travel times from the injector to each producer should be
% equal. When using a skew grid that is not K-orhtogonal, the travel times
% will not be equal and the flow pattern will differ quite a lot from being
% symmetric. In particular, since our discretization method is not
% consistent, the dissymmetry does not decay with increasing grid
% resolution and hence the method does not converge.
mrstModule add incomp coarsegrid msrsb diagnostics streamlines

figure('Position', [440 150 865 351]);

% Rectangular reservoir with a skew grid.
G = cartGrid([41,20],[2,1]);
makeSkew = @(c) c(:,1) + .4*(1-(c(:,1)-1).^2).*(1-c(:,2));
G.nodes.coords(:,1) = 2*makeSkew(G.nodes.coords);
G = computeGeometry(G);
disp(['  Grid: ' num2str(G.cartDims)]);

% Homogeneous reservoir properties
rock = makeRock(G, 100*milli*darcy, .2);
hT   = computeTrans(G, rock);
pv   = sum(poreVolume(G,rock));

% Symmetric well pattern
srcCells = findEnclosingCell(G,[2 .975; .4 .025; 3.6 .025]);
src = addSource([], srcCells, [pv; -.5*pv; -.5*pv]);

% Single-phase fluid
fluid = initSingleFluid('mu', 1*centi*poise,'rho', 1000*kilogram/meter^3);

% Initial state
state  = initState(G,[], 0);

for i=1:2
    if i==1
        p = partitionUI(G, [7 5]);
    else
        p = sampleFromBox(G,reshape(1:35,[7,5]));
    end
    
    % Compute basis functions
    A     = getIncomp1PhMatrix(G, hT, state, fluid);
    CG    = generateCoarseGrid(G, p);
    CG    = coarsenGeometry(CG);
    CG    = storeInteractionRegion(CG);
    basis = getMultiscaleBasis(CG, A, 'type', 'msrsb');
    
    % Compute multiscale solution
    state  = incompTPFA(state, G, hT, fluid, 'src', src);
    ms     = incompMultiscale(state, CG, hT, fluid, basis, 'src', src);
    tof    = computeTimeOfFlight(state, G, rock, 'src', src);
    tofms  = computeTimeOfFlight(ms,    G, rock, 'src', src);
    
    subplot(2,2,i)
    plotCellData(G, ms.pressure, 'EdgeColor', 'none');
    outlineCoarseGrid(G, p, 'Color','k','LineWidth',1);
    hold on
    plot(G.cells.centroids(srcCells,1), G.cells.centroids(srcCells,2),...
        '.','Color',[.9 .9 .9],'MarkerSize',16);
    hold off
    colormap(gca,.85*parula(16).^2+.15);
    
    subplot(2,2,2+i);
    plotCellData(G, tofms, tofms<.2, 'EdgeColor','none'); caxis([0 .2]); box on
    seed = floor(G.cells.num/5)+(1:G.cartDims(1))';
    hf = streamline(pollock(G, state, seed, 'substeps', 1) );
    hb = streamline(pollock(G, state, seed, 'substeps', 1, 'reverse' , true));
    set ([ hf ; hb ], 'Color' , 'k' );
end
%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

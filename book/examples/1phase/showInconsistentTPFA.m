%% Demonstrate lack of convergence for the TPFA scheme
% We consider a homogeneous, rectangular reservoir with a symmetric well
% pattern consisting of one injector and two producers. Because of the
% symmetry, the travel times from the injector to each producer should be
% equal. When using a skew grid that is not K-orhtogonal, the travel times
% will not be equal and the flow pattern will differ quite a lot from being
% symmetric. In particular, since our discretization method is not
% consistent, the dissymmetry does not decay with increasing grid
% resolution and hence the method does not converge.
mrstModule add incomp diagnostics streamlines

figure('Position', [440 450 865 351]);
T = nan(30,2);
disp('Convergence study:');
for i=1:1:30
   % Rectangular reservoir with a skew grid.
   G = cartGrid([i*20+1,i*10],[2,1]);
   makeSkew = @(c) c(:,1) + .4*(1-(c(:,1)-1).^2).*(1-c(:,2));
   G.nodes.coords(:,1) = 2*makeSkew(G.nodes.coords);
   G = computeGeometry(G);
   disp(['  Grid: ' num2str(G.cartDims)]);
   
   % Homogeneous reservoir properties
   rock = makeRock(G, 100*milli*darcy, .2);
   pv   = sum(poreVolume(G,rock));

   % Symmetric well pattern
   srcCells = findEnclosingCell(G,[2 .975; .5 .025; 3.5 .025]);
   src = addSource([], srcCells, [pv; -.5*pv; -.5*pv]);

   % Single-phase fluid
   fluid = initSingleFluid('mu', 1*centi*poise,'rho', 1000*kilogram/meter^3);

   % Solve flow problem
   hT     = computeTrans(G, rock);
   state  = initState(G,[], 0);
   state  = incompTPFA(state, G, hT, fluid, 'src', src);
   tof    = computeTimeOfFlight(state, G, rock, 'src', src);
   T(i,:) = tof(srcCells(2:3))';

   % Plot second solution
   if i==2
      subplot(2,3,1:2);
      plotCellData(G, state.pressure, 'EdgeColor', 'k', 'EdgeAlpha', .05);
      hold on
      plot([.5 2 3.5], [.025 .975 .025],'.','Color',[.9 .9 .9],'MarkerSize',16);
      hold off

      subplot(2,3,4:5);
      plotCellData(G, tof, tof<.2, 'EdgeColor','none'); caxis([0 .2]); box on
      seed = floor(G.cells.num/5)+(1:G.cartDims(1))';
      hf = streamline(pollock(G, state, seed, 'substeps', 1) );
      hb = streamline(pollock(G, state, seed, 'substeps', 1, 'reverse' , true));
      set ([ hf ; hb ], 'Color' , 'k' );
      drawnow;
   end
end
subplot(2,3,[3 6]);
bar(T')
hold on, plot([.5 2.5],[1 1], '--r','LineWidth',2); hold off
set(gca,'XTickLabel', {'L','R'});
axis tight

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

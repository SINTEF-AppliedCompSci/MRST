%% Norne: Two-Phase Incompressible Simulator
% In this example we will solve an incompressible two-phase oil-water
% problem, which consists of an elliptic pressure equation
%
% $$\nabla\cdot v = q, \qquad v=\textbf{--}\lambda K\nabla p,$$
%
% where v is the Darcy velocity (total velocity) and lambda is the
% mobility, which depends on the water saturation S.
%
% The saturation equation (conservation of the water phase) is given as:
%
% $$ \phi \frac{\partial S}{\partial t} +
%     \nabla \cdot (f_w(S) v) = q_w$$
%
% where phi is the rock porosity, f is the Buckley-Leverett fractional flow
% function, and q_w is the water source term.
%
% This is an independent continuation of
% <maltab:edit('incompExampleNorne1ph.m') incompExampleNorne1ph>, in
% which we solved the corresponding single-phase problem for the
% <matlab:exit('showNorne.m') Norne model>, which is a real field from the
% Norwegian Sea.

mrstModule add incomp
linsolve = @mldivide;

%% Setup the model
% How to read and setup the geological model was discussed in detail in the
% <matlab:exit('showNorne.m') showNorne> and
% <maltab:edit('incompExampleNorne1ph.m') incompExampleNorne1ph> tutorials.
% Here, we therefore just exectute the necessarry commands without further
% comments,
grdecl = fullfile(getDatasetPath('norne'), 'NORNE.GRDECL');
grdecl = readGRDECL(grdecl);
grdecl = convertInputUnits(grdecl, getUnitSystem('METRIC'));

G = processGRDECL(grdecl); clear grdecl;
G = computeGeometry(G(1));

%% Set rock and fluid data
% The permeability is lognormal and isotropic within nine distinct layers
% and is generated using our simplified 'geostatistics' function and then
% transformed to lay in the interval 200-2000 mD. For the
% permeability-porosity relationship we use the simple relationship that
% phi~0.25*K^0.11, porosities in the interval 0.25-0.32. For the two-phase
% fluid model, we use values:
%
% * densities: [rho_w, rho_o] = [1000 700] kg/m^3
% * viscosities: [mu_w, mu_o] = [1 5] cP.
gravity off
K          = logNormLayers(G.cartDims, rand([9, 1]), 'sigma', 2);
K          = K(G.cells.indexMap);
K          = 200 + (K-min(K))/(max(K)-min(K))*1800;
rock.perm  = convertFrom(K, milli*darcy);
rock.poro  = 0.25 * (K ./ 200).^0.1;
rock.ntg   = ones(size(K)); 
fluid      = initSimpleFluid('mu' , [   1,   5]*centi*poise     , ...
                             'rho', [1014, 859]*kilogram/meter^3, ...
                             'n'  , [   2,   2]);

%% Visualize the model
% Set various visualization parameters
myview = struct('vw',   [90, 20],    ...  % view angle
                'zoom', 1,        ...  % zoom
                'asp',  [1 1 0.2],  ...  % data aspect ratio
                'wh',   30,         ...  % well height above reservoir
                'cb',   'horiz'     ...  % colorbar location
                );

% Plot the data
clf, title('Log_{10} of x-permeability [mD]');
plotCellData(G,log10(rock.perm),'EdgeColor','k','EdgeAlpha',0.1);
set(gca,'dataasp',myview.asp);
axis off, view(myview.vw); zoom(myview.zoom), colormap(jet), axis tight

cs = [200 400 700 1000 1500 2000];
caxis(log10([min(cs) max(cs)]*milli*darcy));
h = colorbarHist(log10(rock.perm),caxis,'South',100);
set(h, 'XTick', log10(cs*milli*darcy), 'XTickLabel', num2str(cs'));

%% Introduce wells
% The reservoir is produced using a set production wells controlled by
% bottom-hole pressure and rate-controlled injectors. Wells are described
% using a Peacemann model, giving an extra set of equations that need to be
% assembled. For simplicity, all wells are assumed to be vertical and are
% assigned using the logical (i,j) subindex.

% Set vertical injectors, completed in the lowest 12 layers.
nz = G.cartDims(3);
I = [ 9, 26,  8, 25, 35, 10];
J = [14, 14, 35, 35, 68, 75];
R = [ 4,  4,  4,  4,  4,  4]*1000*meter^3/day;
nIW = 1:numel(I); W = [];
for i = 1 : numel(I)
   W = verticalWell(W, G, rock, I(i), J(i), nz-11:nz, 'Type', 'rate', ...
                    'InnerProduct', 'ip_tpf', ...
                    'Val', R(i), 'Radius', 0.1, 'Comp_i', [1, 0], ...
                    'name', ['I$_{', int2str(i), '}$'], ...
                    'RefDepth', 2400*meter);  % Above formation
end

% Set vertical producers, completed in the upper 14 layers
I = [17, 12, 25, 35, 15];
J = [23, 51, 51, 95, 94];
P = [300, 300, 300, 200, 200];
nPW = (1:numel(I))+max(nIW);
for i = 1 : numel(I)
   W = verticalWell(W, G, rock, I(i), J(i), 1:14, 'Type', 'bhp', ...
                    'InnerProduct', 'ip_tpf', ...
                    'Val', 300*barsa(), 'Radius', 0.1, ...
                    'name', ['P$_{', int2str(i), '}$'], ...
                    'Comp_i', [0, 1], 'RefDepth', 2400*meter);
end

% Plot grid outline and the wells
clf
plotGrid(G,'FaceColor','none','EdgeAlpha',0.1);
axis tight off,
view(myview.vw), set(gca,'dataasp',myview.asp), zoom(myview.zoom*.75);
plotWell(G,W,'height',50);
plotGrid(G, vertcat(W(nIW).cells), 'FaceColor', 'b');
plotGrid(G, vertcat(W(nPW).cells), 'FaceColor', 'r');
camdolly(0, .3, 0);


%% Transmissibilities and initial state
% Initialize solution structures and compute transmissibilities from
% input grid, rock properties, and well structure.
trans = computeTrans(G, rock, 'Verbose', true);
rSol  = initState(G, W, 0, [0, 1]);

%% Solve initial pressure
% Solve linear system construced from S and W to obtain solution for flow
% and pressure in the reservoir and the wells.
rSol = incompTPFA(rSol, G, trans, fluid, 'wells', W, ...
                       'LinSolve', linsolve);

clf, title('Initial pressure')
plotCellData(G, convertTo(rSol.pressure(1:G.cells.num), barsa), ...
    'EdgeColor', 'k', 'EdgeAlpha', 0.1);
plotWell(G, W, 'height', myview.wh, 'color', 'k');
axis tight off, view(myview.vw)
set(gca,'dataasp',myview.asp), zoom(myview.zoom*.9)
colorbar(myview.cb), colormap(jet)
camdolly(0,.3,0)

%% Main loop
% In the main loop, we alternate between solving the transport and the flow
% equations. The transport equation is solved using the standard implicit
% single-point upwind scheme with a simple Newton-Raphson nonlinear solver.
T      = 6*year();
dT     = .5*year;
dTplot = dT;
pv     = poreVolume(G,rock);

% Prepare plotting of saturations
clf
plotGrid(G, 'FaceColor', 'none', 'EdgeAlpha', .1);
plotWell(G, W, 'height', myview.wh, 'color', 'k');
view(myview.vw), set(gca,'dataasp',myview.asp), zoom(myview.zoom/1.2);
axis tight off
h=colorbar(myview.cb); colormap(flipud(winter)),
set(h,'Position',[.13 .07 .77 .05]);
camdolly(0,0.2,0)
[hs,ha] = deal([]); caxis([0 1]);

%% Start the main loop
t  = 0;  plotNo = 1;
while t < T,
   rSol = implicitTransport(rSol, G, dT, rock, fluid, 'wells', W);

   % Check for inconsistent saturations
   assert(max(rSol.s(:,1)) < 1+eps && min(rSol.s(:,1)) > -eps);

   % Update solution of pressure equation.
   rSol = incompTPFA(rSol, G, trans, fluid, 'wells', W);

    % Increase time and continue if we do not want to plot saturations
   t = t + dT;
   if ( t < plotNo*dTplot && t <T), continue, end

   % Plot saturation
   delete([hs, ha])
   hs = plotCellData(G, rSol.s(:,1), find(rSol.s(:,1) > 0.01),'EdgeColor','none');
   ha = annotation('textbox', [0 0.93 0.32 0.07], ...
                   'String', ['Water saturation at ', ...
                              num2str(convertTo(t,year)), ' years'],'FontSize',8);
   drawnow
   plotNo = plotNo+1;
end

%%
displayEndOfDemoMessage(mfilename)

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
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

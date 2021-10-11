%% SAIGUP: Solving Two-Phase Flow on a Realistic Corner-Point Model
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
% <maltab:edit('incompExampleSAIGUP1ph.m') incompExampleSAIGUP1ph>, in
% which we solved the corresponding single-phase problem for the
% <matlab:exit('showSAIGUP.m') SAIGUP model>, which is a synthetic, but
% realistic model of a shallow-marine reservoir.

mrstModule add incomp

%% Decide which linear solver to use
%
% We use MATLAB(R)'s built-in <matlab:doc('mldivide') |mldivide|>
% ("backslash") linear solver software to resolve all systems of linear
% equations that arise in the various discretzations.  In many cases, a
% multigrid solver such as Notay's AGMG package may be preferable.  We
% refer the reader to Notay's home page at
% http://homepages.ulb.ac.be/~ynotay/AGMG/ for additional details on this
% package.
linsolve_p = @mldivide;  % Pressure
linsolve_t = @mldivide;  % Transport (implicit)

% If you have AGMG available
% mrstModule add agmg
% linsolve_p = @agmg;  % Pressure
% linsolve_t = @agmg;  % Transport (implicit)

%% Setup the model
% How to read and setup the geological model was discussed in detail in the
% <matlab:exit('showSAIGUP.m') showSAIGUP> and
% <maltab:edit('incompExampleSAIGUP1ph.m') incompExampleSAIGUP1ph>
% tutorials. Here, we therefore just exectute the necessarry commands
% without further comments,
grdecl = fullfile(getDatasetPath('SAIGUP'), 'SAIGUP.GRDECL');
grdecl = readGRDECL(grdecl);
grdecl = convertInputUnits(grdecl, getUnitSystem('METRIC'));

G = processGRDECL  (grdecl);
G = computeGeometry(G);

rock = grdecl2Rock(grdecl, G.cells.indexMap);
is_pos                = rock.perm(:, 3) > 0;
rock.perm(~is_pos, 3) = min(rock.perm(is_pos, 3));

%%
% Plot the logarithm of the permeability in the z-direction. The model has
% a large number of very small permeabilities (shown in dark blue), and in
% the histogram we have only included permeability values that are larger
% than 1 micro darcy.
myview = struct('vw',   [-110,18],  ...  % view angle
                'zoom', 1.0,        ...  % zoom
                'asp',  [15 15 2],  ...  % data aspect ratio
                'wh',   50,         ...  % well height above reservoir
                'cb',   'horiz'     ...  % colorbar location
                );
K   = rock.perm(:,3);  idx = K>1e-3*milli*darcy;
clf,
plotCellData(G,log10(K),'EdgeColor','k','EdgeAlpha',0.1);
set(gca,'dataasp',myview.asp);
axis off, view(myview.vw); zoom(myview.zoom), colormap(jet)

cs = [0.001 0.01 0.1 1 10 100 1000];
caxis(log10([min(cs) max(cs)]*milli*darcy));
h = colorbarHist(log10(K(idx)),caxis,'South',100);
set(h, 'XTick', log10(cs*milli*darcy), 'XTickLabel', num2str(cs'));
camdolly(.2, 0, 0);
 
%% Set fluid data
% For the two-phase fluid model, we use values:
%
% * densities: [rho_w, rho_o] = [1000 700] kg/m^3
% * viscosities: [mu_w, mu_o] = [1 5] cP.
gravity off
fluid      = initSimpleFluid('mu' , [   1,   5]*centi*poise     , ...
                             'rho', [1000, 700]*kilogram/meter^3, ...
                             'n'  , [   2,   2]);

%% Introduce wells
% The reservoir is produced using a set of production wells controlled by
% bottom-hole pressure and rate-controlled injectors. Wells are described
% using a Peacemann model, giving an extra set of equations that need to be
% assembled. For simplicity, all wells are assumed to be vertical and are
% assigned using the logical (i,j) subindex.

% Set vertical injectors, completed in the lowest 12 layers.
nz = G.cartDims(3);
I = [ 9,  8, 25, 25, 10];
J = [14, 35, 35, 95, 75];
R = [ 4,  4,  4,  4,  4,  4]*125*meter^3/day;
nIW = 1:numel(I); W = [];
for i = 1 : numel(I),
   W = verticalWell(W, G, rock, I(i), J(i), nz-11:nz, 'Type', 'rate', ...
                    'InnerProduct', 'ip_tpf', ...
                    'Val', R(i), 'Radius', 0.1, 'Comp_i', [1, 0], ...
                    'name', ['I$_{', int2str(i), '}$']);
end

% Set vertical producers, completed in the upper 14 layers
I = [17, 12, 25, 33, 7];
J = [23, 51, 51, 95, 94];
nPW = (1:numel(I))+max(nIW);
for i = 1 : numel(I),
   W = verticalWell(W, G, rock, I(i), J(i), 1:14, 'Type', 'bhp', ...
                    'InnerProduct', 'ip_tpf', ...
                    'Val', 300*barsa(), 'Radius', 0.1, ...
                    'name', ['P$_{', int2str(i), '}$'], 'Comp_i', [0, 1]);
end

% Plot grid outline and the wells
clf
plotGrid(G,'FaceColor','none','EdgeAlpha',0.1);
axis tight off,
view(myview.vw), set(gca,'dataasp',myview.asp), zoom(myview.zoom);
plotWell(G,W,'height',50);
plotGrid(G, vertcat(W(nIW).cells), 'FaceColor', 'b');
plotGrid(G, vertcat(W(nPW).cells), 'FaceColor', 'r');

%% Compute transmissibilities and init reservoir
trans = computeTrans(G, rock, 'Verbose', true);
rSol  = initState(G, W, 0, [0, 1]);

%% Construct pressure and transport solvers
solve_press  = @(x) incompTPFA(x, G, trans, fluid, 'wells', W, ...
                                    'LinSolve', linsolve_p);
solve_transp = @(x, dt) ...
   implicitTransport(x, G, dt, rock, fluid, ...
                     'wells', W, 'LinSolve', linsolve_t);

%% Solve initial pressure
% Solve linear system construced from S and W to obtain solution for flow
% and pressure in the reservoir and the wells.
rSol = solve_press(rSol);

clf, title('Initial pressure')
plotCellData(G, convertTo(rSol.pressure(1:G.cells.num), barsa), ...
    'EdgeColor', 'k', 'EdgeAlpha', 0.1);
plotWell(G, W, 'height', myview.wh, 'color', 'k');
axis tight off, view(myview.vw)
set(gca,'dataasp',myview.asp), zoom(myview.zoom), camdolly(.1,.15,0)
caxis([300 400]), colorbar(myview.cb), colormap(jet)

%% Main loop
% In the main loop, we alternate between solving the transport and the flow
% equations. The transport equation is solved using the standard implicit
% single-point upwind scheme with a simple Newton-Raphson nonlinear solver.
T      = 12*year();
dT     = year();
dTplot = year();
pv     = poreVolume(G,rock);

% Prepare plotting of saturations
clf
plotGrid(G, 'FaceColor', 'none', 'EdgeAlpha', .1);
plotWell(G, W, 'height', myview.wh, 'color', 'k');
view(myview.vw), set(gca,'dataasp',myview.asp), zoom(myview.zoom);
axis tight off
h=colorbar(myview.cb); colormap(flipud(winter)),
set(h,'Position',[.13 .07 .77 .05]);
camdolly(.05,0,0)
[hs,ha] = deal([]); caxis([0 1]);

% Start the main loop
t  = 0;  plotNo = 1;
while t < T,
   rSol = solve_transp(rSol, dT);

   % Check for inconsistent saturations
   assert(max(rSol.s(:,1)) < 1+eps && min(rSol.s(:,1)) > -eps);

   % Update solution of pressure equation.
   rSol = solve_press(rSol);

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

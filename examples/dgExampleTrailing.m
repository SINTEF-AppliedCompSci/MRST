%% Resolution of trailing linear waves
% In this example, we consider an incompressible Buckley-Leverett
% displacement with an approximate, piecewise linear flux function as an
% idealized model of trailing linear waves seen in multicomponent and
% compositional simulation models. The example shows that dG(1) gives
% similar resolution of the leading displacement wave as an explicit SPU
% scheme and better resolution of the trailing linear waves.

%% Add modules
mrstModule add dg ad-core ad-props ad-blackoil sequential

%% Set up problem
% Construct grid, compute geometry and cell dimensions, and set
% petrophysical properties 
n    = 40; m=25/40*n;
G    = computeGeometry(cartGrid([n,1],[1 1]));
G    = computeCellDimensions(G);
rock = makeRock(G, 1, 1);

% We consider Bucley-Leverett-type displacement with quadratic relperms
fluid = initSimpleADIFluid('phases', 'WO' , 'n', [2,3], ... 
                           'mu', [1,1],  'rho', [1,1]);
sm = linspace(0,1,6);
fluid.f_w = @(s) interpTable(sm, sm.^2./(sm.^2 + (1-sm).^3), s);

figure
fm = fluid.f_w(sm);
sc  = sm([end:-1:4 1]);
fc  = fm([end:-1:4 1]);
dfc = diff(fc)./diff(sc); 
plot(sm,fm,'-o'); hold on,plot(sc,fc,'--','LineWidth',1); hold off
legend('f(S)','f^c(S)');

% The base model is a generic black-oil model with oil and water
model  = GenericBlackOilModel(G, rock, fluid, 'gas', false);

% Initial state: filled with oil and unit volumetric flow rate
state0 = initResSol(G, 1, [0,1]);
state0.flux(1:G.cells.num+1) = 1;

% Boundary conditions
bc = fluxside([], G, 'left' ,  1, 'sat', [1,0]); % Inflow
bc = fluxside(bc, G, 'right', -1, 'sat', [1,0]); % Outflow

% Define schedule: unit time steps, simulate almost to breakthrough
schedule = simpleSchedule(ones(m,1)./n, 'bc', bc);

%% Finite volume simulation
figure, 
plot([0 dfc(rldecode((1:numel(dfc))',2))*.625 1],sc(rldecode((1:numel(sc))',2)));
hold all

% Use standard discretization in MRST: single-point upwind (SPU)
tmodel = TransportModel(model);
tmodel = tmodel.validateModel();
tmodel.parentModel.FlowDiscretization.ComponentPhaseFlux = ComponentPhaseFluxFractionalFlowSimple(model);
[~, stSPUi, repSPUi] = simulateScheduleAD(state0, tmodel, schedule);
plot(G.cells.centroids(:,1),stSPUi{end}.s(:,1),'.','MarkerSize',20);

%% dG implicit
tmodelDGi = TransportModelDG(model, 'degree', [1,0]);
tmodelDGi.storeUnlimited = true;
tmodelDGi = tmodelDGi.validateModel();
tmodelDGi.parentModel.FlowDiscretization.ComponentPhaseFlux = ComponentPhaseFluxFractionalFlowSimple(model);
tmodelDGi.parentModel.FlowDiscretization.ComponentPhaseVelocity = ComponentPhaseVelocityFractionalFlowSimpleDG(model);
[~, stDGi, repDGi] = simulateScheduleAD(state0, tmodelDGi, schedule);
plot(G.cells.centroids(:,1),stDGi{end}.s(:,1),'.','MarkerSize',20);

%% SPU explicit
tmodel3 = tmodel.validateModel;
flux = tmodel3.parentModel.FlowDiscretization;
fb   = ExplicitFlowStateBuilder();
%fb    = AdaptiveImplicitFlowStateBuilder();
flux = flux.setFlowStateBuilder(fb);
tmodel3.parentModel.FlowDiscretization = flux;
[~, stSPUe, repSPUe] = simulateScheduleAD(state0, tmodel3, schedule);
plot(G.cells.centroids(:,1),stSPUe{end}.s(:,1),'.','MarkerSize',20);
 
%% dG explicit
% Not supported yet!

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

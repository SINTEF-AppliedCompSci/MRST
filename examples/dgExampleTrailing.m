%% Bucley-Leverett displacement
% In this example, we consider an incompressible Buckley-Leverett
% displacement in a horizontal 1D channel aligned with the x-axis, and
% compare higher-order discontinuous Galerkin methods with the standard
% finite-volume discretization.

%% Add modules
mrstModule add dg ad-core ad-props ad-blackoil blackoil-sequential

%% Set up problem
% Construct grid, compute geometry and cell dimensions, and set
% petrophysical properties 
n    = 80; m=25/40*n;
G    = computeGeometry(cartGrid([n,1],[1 1]));
G    = computeCellDimensions(G);
rock = makeRock(G, 1, 1);

% We consider Bucley-Leverett-type displacement with quadratic relperms
fluid = initSimpleADIFluid('phases', 'WO' , 'n', [2,2], ... 
                           'mu', [1,1],  'rho', [1,1]);
sm = linspace(0,1,6);
% sm = [0:0.2:0.6 linspace(0.8,1,21)];
fluid.krO = @(s) interpTable(sm, sm.^3, s);
fluid.krW = @(s) interpTable(sm, sm.^2, s);
fluid.f_w = @(s) interpTable(sm, sm.^2./(sm.^2 + (1-sm).^3), s);
fm = fluid.f_w(sm);
sc  = sm([end:-1:4 1]);
fc  = fm([end:-1:4 1]);
dfc = diff(fc)./diff(sc); 
plot(sm,fm,'-o'); hold on,plot(sc,fc,'--','LineWidth',1); hold off

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
tmodel.parentModel.FluxDiscretization.ComponentPhaseFlux = ComponentPhaseFluxFractionalFlowSimple(model);
[~, stSPUi, repSPUi] = simulateScheduleAD(state0, tmodel, schedule);
plot(G.cells.centroids(:,1),stSPUi{end}.s(:,1),'.','MarkerSize',20);

%% dG implicit
tmodelDGi = TransportModelDG(model, 'degree', [1,0]);
tmodelDGi = tmodelDGi.validateModel();
tmodelDGi.parentModel.FluxDiscretization.ComponentPhaseFlux = ComponentPhaseFluxFractionalFlowSimple(model);
tmodelDGi.parentModel.FluxDiscretization.ComponentPhaseVelocity = ComponentPhaseVelocityFractionalFlowSimpleDG(model);
[~, stDGi, repDGi] = simulateScheduleAD(state0, tmodelDGi, schedule);
plot(G.cells.centroids(:,1),stDGi{end}.s(:,1),'.','MarkerSize',20);

%% SPU explicit
tmodel3 = tmodel.validateModel;
flux = tmodel3.parentModel.FluxDiscretization;
fb   = ExplicitFlowStateBuilder();
%fb    = AdaptiveImplicitFlowStateBuilder();
flux = flux.setFlowStateBuilder(fb);
tmodel3.parentModel.FluxDiscretization = flux;
[~, stSPUe, repSPUe] = simulateScheduleAD(state0, tmodel3, schedule);
plot(G.cells.centroids(:,1),stSPUe{end}.s(:,1),'.','MarkerSize',20);

%% dG explicit
tmodelDGe = tmodelDGi;
flux      = tmodelDGi.parentModel.FluxDiscretization;
fb        = ExplicitFlowStateBuilderDG();
flux      = flux.setFlowStateBuilder(fb);
tmodelDGe.parentModel.FluxDiscretization = flux;
[~, stDGe, repDGe] = simulateScheduleAD(state0, tmodelDGe, schedule);
plot(G.cells.centroids(:,1),stSPUe{end}.s(:,1),'.','MarkerSize',20);

%% WENO explicit
tmodel4 = tmodel3;
tmodel4.parentModel.operators.faceUpstr = @(flag, x) weno.faceUpstream(model, [], flag, x);
[~, stWENOe, repWENOe] = simulateScheduleAD(state0, tmodel4, schedule);
plot(G.cells.centroids(:,1),stWENOe{end}.s(:,1),'.','MarkerSize',20);

hold off

%% First-order dG simulation
tmodelDG1 = TransportModelDG(model, 'degree', [1,0]);
%tmodelDG1 = tmodelDG1.validateModel();
%tmodelDG1.parentModel.FluxDiscretization.ComponentPhaseVelocity = ComponentPhaseFluxFractionalFlowSimple(model);
[~, stDG1, repDG1] = simulateScheduleAD(state0, tmodelDG1, schedule);

%% Higher-order dG simulations
tmodelDG2 = TransportModelDG(model, 'degree', [2,0]);
[~, stDG2, repDG2] = simulateScheduleAD(state0, tmodelDG2, schedule);
                                                    

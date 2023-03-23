%% Two-phase gas injection in a water-saturated dual-porosity reservoir. 
% Gas flows in the fracture system, and transfer to the matrix happens via
% gravity drainage.
clear;
clc;
close all;

mrstModule add ad-blackoil ad-core ad-props dual-porosity

%% Set up grid
G = cartGrid([60, 1, 10], [100, 1, 100]*meter);
G = computeGeometry(G);

%% Set up rock properties
rock_matrix = makeRock(G, 10*milli*darcy, .1);
rock_fracture = makeRock(G, 100*milli*darcy, .01);

%% Set up fluid
fluid_matrix = initSimpleADIFluid('phases', 'WO', 'c', [1e-12;1e-12],...
                                  'rho',[1000,1]);
fluid_fracture = initSimpleADIFluid('phases', 'WO', 'c', [1e-12;1e-12],...
                                  'rho',[1000,1]);

%% Set the DP model. Here, a two-phase model is used. Rock and fluid
% are sent in the constructor as a cell array {fracture,matrix}.  
% We use and OilWater model, where oil plays the role of the gas.
gravity on;
model = TwoPhaseOilWaterModelDP(G, {rock_fracture,rock_matrix},...
                                      {fluid_fracture,fluid_matrix});                                  
                     
%% Setting transfer function. This step is important to ensure that fluids
% will be transferred from fracture to matrix (and vice-versa). There are 
% several options of shape factor (see folder
% transfer_models/shape_factor_models/) that could be set in this
% constructor. The second argument is the matrix block size.
model.transfer_model_object = EclipseTwoPhaseTransferFunction('CoatsShapeFactor',[1,1,10]);
                       
%% Initializing state                       
state = initResSol(G, 0*psia, [1,0]);
state.wellSol = initWellSolAD([], model, state);
state.pressure_matrix = state.pressure;
state.sm = state.s;

%% Boundary conditions
bc = pside([], G, 'xmax', 0*psia, 'sat', [0,1]);
bc = pside(bc, G, 'xmin', 500*psia, 'sat', [0,1]);

%% Solver
solver = NonLinearSolver();

%% Validating model
model = model.validateModel();

%% Figure
fig1 = figure('Position',[100,100,1200,600]);
fig1.Color = 'w';
colormap('jet');

%% Time loop
dt = 5*day;
tmax = 100*dt;
t = 0;
while t<=tmax
    
    disp(['Time = ',num2str(t/day), ' days'])
    state = solver.solveTimestep(state, dt, model, 'bc', bc);

    set(0, 'CurrentFigure', fig1);

    subplot(2,1,1);
    p = plotCellData(G,state.s(:,1));
    p.EdgeAlpha = 0;
    view(0,0);
    colorbar;
    caxis([0,1]);
    set(gca,'FontSize',16);
    xlabel('x')
    zlabel('z')
    title('Fractures');
    
    subplot(2,1,2);
    p = plotCellData(G,state.sm(:,1));
    p.EdgeAlpha = 0;
    view(0,0);
    colorbar;
    caxis([0,1]);
    set(gca,'FontSize',16);
    xlabel('x')
    zlabel('z')
    title('Matrix');

    drawnow;

    t = t+dt;
end

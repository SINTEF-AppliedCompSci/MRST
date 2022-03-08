%% Single-phase depletion of a dual-porosity reservoir. Fracture fluid is
% produced before the matrix. Note that the pressure drop is faster in the
% fractures, while the matrix drains slowly.
clear;
clc;
close all;

mrstModule add ad-blackoil ad-core ad-props dual-porosity

%% Set up grid
G = cartGrid([100, 1, 1], [60, 1, 1]*meter);
G = computeGeometry(G);

%% Set up rock properties
rock_matrix = makeRock(G, 1*darcy, .3);
rock_fracture = makeRock(G, 1000*darcy, .01);

%% Set up fluid
fluid_matrix = initSimpleADIFluid('phases', 'W', 'c', 1e-06);
fluid_fracture = initSimpleADIFluid('phases', 'W', 'c', 1e-06);

%% Set the DP model. Here, a single-phase model is used. Rock and fluid
% are sent in the constructor as a cell array {fracture,matrix}
model = WaterModelDP(G, {rock_fracture,rock_matrix},...
                        {fluid_fracture,fluid_matrix});
                     
%% Setting transfer function. This step is important to ensure that fluids
% will be transferred from fracture to matrix (and vice-versa). There are 
% several options of shape factor (see folder
% transfer_models/shape_factor_models/) that could be set in this
% constructor. The second argument is the matrix block size.
model.transfer_model_object = KazemiSinglePhaseTransferFunction('KazemiShapeFactor',...
                                                                [10,10,10]);
                       
%% Initializing state                       
state = initResSol(G, 1000*psia);
state.wellSol = initWellSolAD([], model, state);
state.pressure_matrix = state.pressure;

%% Boundary conditions
bc = pside([], G, 'xmin', 0*psia, 'sat', 1);

%% Solver
solver = NonLinearSolver();

%% Validating model
model = model.validateModel();

%% Figure
fig1 = figure('Position',[100,100,1200,600]);
fig1.Color = 'w';
colormap('jet');

%% Time loop
dt = 0.1*day;
tmax = 100*dt;
t = 0;
while t<=tmax
    
    disp(['Time = ',num2str(t/day), ' days'])
    state = solver.solveTimestep(state, dt, model, 'bc', bc);

    set(0, 'CurrentFigure', fig1)

    subplot(2,2,1);
    p = plotCellData(G,state.pressure_matrix/psia);
    p.EdgeAlpha = 0;
    colorbar;
    caxis([0,1000]);
    set(gca,'FontSize',16);
    xlabel('x')
    ylabel('y')

    subplot(2,2,3);
    p = plotCellData(G,state.pressure/psia);
    p.EdgeAlpha = 0;
    colorbar;
    caxis([0,1000]);
    set(gca,'FontSize',16);
    xlabel('x')
    ylabel('y')

    subplot(2,2,2);
    plot(G.cells.centroids(:,1),state.pressure_matrix/psia,'LineWidth',1.5,'Color','r');
    set(gca,'FontSize',16);
    xlabel('x')
    ylim([0,1000])
    ylabel('Pm [psia]')

    subplot(2,2,4);
    plot(G.cells.centroids(:,1),state.pressure/psia,'LineWidth',1.5);
    set(gca,'FontSize',16);
    xlabel('x')
    ylim([0,1000])
    ylabel('Pf [psia]')

    drawnow;

    t = t+dt;
end

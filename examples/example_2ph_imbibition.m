%% Two-phase water injection in an oil-saturated dual-porosity reservoir. 
% Water flows quickly in the fracture system, while transfer to the matrix
% happens via spontaneous imbibition.
clear;
clc;
close all;

mrstModule add ad-blackoil ad-core ad-props dual-porosity

%% Set up grid
G = cartGrid([60, 1, 1], [60, 1, 1]*meter);
G = computeGeometry(G);

%% Set up rock properties
rock_matrix = makeRock(G, 10*milli*darcy, .1);
rock_fracture = makeRock(G, 100*milli*darcy, .01);

%% Set up fluid
fluid_matrix = initSimpleADIFluid('phases', 'WO', 'c', [1e-12;1e-12]);
fluid_fracture = initSimpleADIFluid('phases', 'WO', 'c', [1e-12;1e-12]);
Pe = 100*psia;
fluid_matrix.pcOW = @(sw)Pe;

%% Set the DP model. Here, a two-phase model is used. Rock and fluid
% are sent in the constructor as a cell array {fracture,matrix}
model = TwoPhaseOilWaterModelDP(G, {rock_fracture,rock_matrix},...
                                   {fluid_fracture,fluid_matrix});
                     
%% Setting transfer function. This step is important to ensure that fluids
% will be transferred from fracture to matrix (and vice-versa). There are 
% several options of shape factor (see folder
% transfer_models/shape_factor_models/) that could be set in this
% constructor. The second argument is the matrix block size. Another
% possible transfer function to be used in this model would be:
%       model.transfer_model_object = EclipseTransferFunction();
model.transfer_model_object = KazemiTwoPhaseTransferFunction('KazemiShapeFactor',...
                                                              [5,5,5]);
                       
%% Initializing state                       
state = initResSol(G, 0*psia, [0,1]);
state.wellSol = initWellSolAD([], model, state);
state.pressure_matrix = state.pressure;
state.sm = state.s;

%% Boundary conditions
bc = pside([], G, 'xmax', 0*psia, 'sat', [1,0]);
bc = pside(bc, G, 'xmin', 1000*psia, 'sat', [1,0]);

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
    
    set(0, 'CurrentFigure', fig1)

    subplot(2,2,1);
    p = plotCellData(G,state.s(:,1));
    p.EdgeAlpha = 0;
    colorbar;
    caxis([0,1]);
    set(gca,'FontSize',16);
    xlabel('x')
    ylabel('y')

    subplot(2,2,3);
    p = plotCellData(G,state.sm(:,1));
    p.EdgeAlpha = 0;
    colorbar;
    caxis([0,1]);
    set(gca,'FontSize',16);
    xlabel('x')
    ylabel('y')

    set(0, 'CurrentFigure', fig1)

    subplot(2,2,2);
    plot(G.cells.centroids(:,1),state.s(:,1),'LineWidth',1.5,'Color','r');
    set(gca,'FontSize',16);
    xlabel('x')
    ylim([0,1])
    ylabel('Swf [-]')

    subplot(2,2,4);
    plot(G.cells.centroids(:,1),state.sm(:,1),'LineWidth',1.5);
    set(gca,'FontSize',16);
    xlabel('x')
    ylim([0,1])
    ylabel('Swm [-]')

    drawnow;

    t = t+dt;
end


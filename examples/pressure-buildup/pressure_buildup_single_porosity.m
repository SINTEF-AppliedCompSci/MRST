%% Single phase pressure buildup during fluid injection 
clear;
clc;
close all;

mrstModule add ad-blackoil ad-core ad-props mrst-gui

%% Set up grid
nx = 50;
ny = 50;
nz = 10;
x_size = 10000; dx = x_size/nx;
y_size = 10000; dy = y_size/ny;
z_size = 100; dz = z_size/nz;
G = cartGrid([nx, ny, nz], [x_size, y_size, z_size]);
G = computeGeometry(G);

%% Set up rock properties
rock = makeRock(G, 1*milli*darcy, .2);

%% Set up fluid
rhow = 1000;
rhoco2 = 600;
fluid = initSimpleADIFluid('phases', 'WO',...
                           'c', [1e-10,0],...
                           'rho', [rhow,rhoco2],...
                           'mu', [1,0.04]*centi*poise);

%% Wells
qinj = 0.016*mega*1e+3/rhoco2/year;
W = [];
W = verticalWell(W, G, rock,  nx/2,   ny/2, (1:10),     ...
               'Type', 'rate', 'Val', qinj, ...
               'Name', 'Inj. Well', 'compi' ,[1,0]);
           
%% Set up model
gravity reset on;
model = TwoPhaseOilWaterModel(G, rock, fluid);

%% Calling validate model
model = model.validateModel();

%% Initializing state 
g = 9.81;
ptop = 20*mega*Pascal;
pinit = ptop + rhow*g*G.cells.centroids(:,3);
state = initResSol(G, pinit, [1,0]);
state.wellSol = initWellSolAD(W, model, state);
   
%% Solver
solver = NonLinearSolver();

%% Figure
fig1 = figure('Position',[100,100,800,600]);
fig1.Color = 'w';
colormap('jet');

%% Plot cells
plot_cells = find(G.cells.centroids(:,2) >= y_size/2);
plot_line = find(G.cells.centroids(:,2) == y_size/2-dy/2 &...
                 G.cells.centroids(:,3) == dz/2);
             
%% Fracture pressure
pfrac = 30*mega*Pascal;

%% Time loop
dt = 1*year;
tmax = 30*dt;
t = 0;
while t<=tmax
    
    disp(['Time = ',num2str(t/year), ' years'])
    state = solver.solveTimestep(state, dt, model, 'W', W);
    
    figure(fig1)
    subplot(2,1,1)
    p = plotCellData(G,state.pressure/(mega*Pascal),plot_cells);
    hold on;
    plotWell(G,W);
    p.EdgeAlpha = 0;
    colorbar;
    caxis([20,35]);
    zlim([0,100]);
    view(-7,50);
    set(gca,'FontSize',20);
    xlabel('x')
    ylabel('y')
    drawnow;
    
    subplot(2,1,2)
    plot(G.cells.centroids(plot_line,1),state.pressure(plot_line)/mega,...
         'LineWidth',1.5,'LineStyle','-','Color','b');
    hold on;
    plot(G.cells.centroids(plot_line,1),pfrac*ones(length(plot_line),1)/mega,...
         'LineWidth',1.5,'LineStyle','--','Color','r');
    text(800,32,'P_{frac}','Color','r','FontSize',16);
    ylim([20,35]); 
    hold off;
    set(gca,'FontSize',20);
    xlabel('x')
    ylabel('Pressure')
    drawnow;
    
    
    t = t+dt;
    
end


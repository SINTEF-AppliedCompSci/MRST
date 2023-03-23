%% Validation of the dfm model against Karimi Fard 2004
clear;
clc;
close all;

%% Modules
mrstModule add ad-blackoil ad-core ad-props mrst-gui

%% Set up grid
nx = 40;
ny = 40;
x_size = 1;
y_size = 1;
G = cartGrid([nx, ny], [x_size, y_size]);
G = computeGeometry(G);

%% Defining fractures
[xf, yf] = deal( G.faces.centroids(:,1),G.faces.centroids(:,2));

frac1 = find(abs(yf-0.2) <= 1e-5 & xf <= 0.6);
frac2 = find(abs(xf-0.3) <= 1e-5 & yf <= 0.4);
frac3 = find(abs(xf-0.7) <= 1e-5 & yf <= 0.7 & yf >= 0.3);

frac_edges = [frac1;frac2;frac3];

%% Set up rock properties
rock_matrix = struct();
rock_matrix.poro = 0.2*ones(G.cells.num,1);
rock_matrix.perm = 1*milli*darcy*ones(G.cells.num,1);

rock_fracture = struct();
rock_fracture.poro = ones(length(frac_edges),1);
rock_fracture.aperture = 0.1*milli*meter * ones(length(frac_edges),1);
rock_fracture.perm = rock_fracture.aperture.^2/12;

%% Set up fluid
fluid_matrix = initSimpleADIFluid('phases', 'WO', 'mu', [1,0.45]*centi*poise);
fluid_fracture = fluid_matrix; 


%% Wells
inj_rate = 0.01*sum(poreVolume(G,rock_matrix))/day;
W = [];
W = verticalWell(W, G, rock_matrix,  1,   1, (1:1),     ...
               'Type', 'rate', 'Val', inj_rate, ...
               'Radius', 0.001,...
               'Name', 'I1', 'compi' ,[1, 0]);  
W = verticalWell(W, G, rock_matrix,  nx,   ny, (1:1),     ...
               'Type', 'bhp', 'Val', 0, ...
               'Radius', 0.001,...
               'Name', 'P1', 'compi' ,[0, 1]);
           
%% Set up model
model = TwoPhaseOilWaterModel(G, rock_matrix, fluid_matrix);

%% Calling validate model
model = model.validateModel();
                                                     
%% Fractured domain manager
manager = FracturedDomainManager();
model = manager.addFracturedDomain(model,...
                                     'dfm',...
                                     frac_edges,...
                                     rock_fracture,...
                                     fluid_fracture);

%% Initializing state (important to use model.G)
state = initResSol(model.G, 0, [0, 1]);
state.wellSol = initWellSolAD(W, model, state);

%% Solver
solver = NonLinearSolver();

%% Figure
fig1 = figure('Position',[100,100,700,600]);
fig1.Color = 'w';
colormap('jet');

%% Time loop
dt = 1*day;
tmax = 200*dt;
t = 0;
pvi = [0.0];
qO = [0.0]; 

while t<=tmax
    
    if abs(pvi(end)-0.0)<1.0e-4 | abs(pvi(end)-0.1)<1.0e-4 |...
            abs(pvi(end)-0.3)<1.0e-4 | abs(pvi(end)-0.5)<1.0e-4
        figure(fig1)
        p = plotCellData(G,state.s(1:nx*ny,1));
        p.EdgeAlpha = 0;
        colorbar;
        xlim([0,x_size]);
        ylim([0,y_size]);
        caxis([0,1]);
        set(gca,'FontSize',20);
        xlabel('x')
        ylabel('y')
        title(['Water saturation ',num2str(t/day), ' days'])
        drawnow;
    end

    disp(['Time = ',num2str((t+dt)/day), ' days'])
    state = solver.solveTimestep(state, dt, model, 'W', W);

    t = t+dt;
    pvi = [pvi inj_rate*t/sum(poreVolume(G,rock_matrix))];
    qO = [qO qO(end)+abs(state.wellSol(2).qOs)*dt/sum(poreVolume(G,rock_matrix))];
    
end

%% Compare oil prodction after 2 pore volume injected
load('OilRecoveryKarimiFard2004.mat');
figure();
hold on;
grid on; 
plot(pvi,qO, 'DisplayName', 'Present work','LineWidth',3,...
    'Color',[0.85 0.33 0.1]);
plot(PVI,FOR, 'DisplayName', 'from Karimi-Fard et al. 2004',...
    'MarkerFaceColor',[0 0 0], 'MarkerEdgeColor',[0 0 0],...
    'Marker','o', 'LineWidth',3, 'LineStyle','none');
xlabel('PV Water Injected');
ylabel('PV Oil Produced');
xlim([0,2]);
ylim([0,1]);
set(gca,'FontSize',20);
legend1 = legend(gca,'show');
set(legend1,'Location','southeast');
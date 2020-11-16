%% Combinaion of dual porosity and DFM models
clear;
clc;
close all;

mrstModule add ad-blackoil ad-core ad-props mrst-gui

%% Set up grid
nx = 50;
ny = 50;
x_size = 1000;
y_size = 1000;
G = cartGrid([nx, ny], [x_size, y_size]);
G = computeGeometry(G);

%% Dual-porosity cells, fault core cells and DFM edges
[xc, xf, yf] = deal(G.cells.centroids(:,1), G.faces.centroids(:,1),...
                    G.faces.centroids(:,2));
dual_porosity_cells = find(((xc >= 400) & (xc <= 500)) | ((xc >= 560) &... 
                                                          (xc <= 650)));
fault_core_cells = find(xc>=500 & xc<=550);    
dfm_edges = find((xf >= 490) & (xf <= 570) & (yf < y_size) & (yf > 0));
             

%% Set up rock properties
base_rock = makeRock(G, 1*milli*darcy, .1);
base_rock.perm(fault_core_cells) = 0;
                
                
base_rock.perm(dual_porosity_cells,:) = 100*milli*darcy;                
base_rock.poro(dual_porosity_cells,:) = 0.1;                

dual_porosity_matrix_rock = struct('perm',1*milli*darcy*ones(length(dual_porosity_cells),1),...
                                   'poro',.1*ones(length(dual_porosity_cells),1));

ap = 0.1*milli*meter;
kf = ap^2/12;
dfm_fracture_rock = struct('perm',kf*ones(length(dfm_edges),1),...
                           'poro',1*ones(length(dfm_edges),1),...
                           'aperture',ap*ones(length(dfm_edges),1));

%% Set up fluid
base_fluid = initSimpleADIFluid('phases', 'WO',...
                                'c', [1e-12,1e-12]/psia,...
                                'n', [1,1],...
                                'mu',[1,1.5]*centi*poise);
dual_porosity_matrix_fluid = base_fluid;
dfm_fracture_fluid = base_fluid;

%% Setting capillary pressure in the matrix of the dual porosity domain
dual_porosity_matrix_fluid.pcOW = @(sw)100*kilo*Pascal*sw.^-0.5;

%% Wells
W = [];
W = verticalWell(W, G, base_rock,  nx,   ny, (1:1), ...
               'Type', 'bhp', 'Val', 500*psia, ...
               'Name', 'P1', 'compi' ,[1,0]);    
W = verticalWell(W, G, base_rock,  1,   1, (1:1), ...
               'Type', 'rate', 'Val', 1e-4, ...
               'Name', 'P1', 'compi' ,[1,0]);     
           
%% Set up model
gravity reset off;
model = TwoPhaseOilWaterModel(G, base_rock, base_fluid);
model = model.validateModel();
                                                     
%% Fractured domain manager
manager = FracturedDomainManager();
transfer_model = KazemiMultiphaseTransferFunction(KazemiShapeFactor([10,10,10]));

model = manager.addFracturedDomain(model,...
                                   'multi_continuum',...
                                   dual_porosity_cells,...
                                   dual_porosity_matrix_rock,...
                                   dual_porosity_matrix_fluid,...
                                   'transfer_model',...
                                   transfer_model);     

model = manager.addFracturedDomain(model,...
                                   'dfm',...
                                   dfm_edges,...
                                   dfm_fracture_rock,...
                                   dfm_fracture_fluid);       
                                                             

%% Matrix cells of dual porosity region for plotting purposes
vcells = model.G.FracturedDomains.domains{1}.virtual_cells;

%% Initializing state 
eps = 1e-6;
state = initResSol(model.G, 1000*psia, [eps,1-eps]);
state.wellSol = initWellSolAD(W, model, state);
   
%% Solver
solver = NonLinearSolver();

%% Figure
fig1 = figure('Position',[100,100,600,600]);
fig1.Color = 'w';
fig2 = figure('Position',[700,100,600,600]);
fig2.Color = 'w';
 

%% Time loop
dt = 10*day;
tmax = 6000*day;
t = 0;
while t<=tmax
    
    if t/day==2500
        disp('save next print')
    end
    
    disp(['Time = ',num2str(t/day), ' days'])
    state = solver.solveTimestep(state, dt, model, 'W', W);
    
    figure(fig1)
    colormap(flipud(jet));
    p = plotCellData(G,state.s(1:nx*ny,1,1));
    p.EdgeAlpha = 0.3;
    colorbar;
    caxis([0,1]);
    set(gca,'FontSize',16);
    xlabel('x')
    ylabel('y')
    hold on;
%     title('Flowing Domain');
    
    figure(fig2)
    colormap(flipud(jet));
    p = plotCellData(G,state.s(vcells,1),dual_porosity_cells);
    p.EdgeAlpha = 0.3;
    colorbar;
    caxis([0,1]);
    xlim([0,x_size]);
    ylim([0,y_size]);
    set(gca,'FontSize',16);
    xlabel('x')
    ylabel('y')
%     title('Stagnant Domain');
    drawnow;
    
    t = t+dt;
    
end


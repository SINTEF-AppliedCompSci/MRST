%% Definition of multirate model with five different constant transfer rates
clear;
clc;
close all;

mrstModule add ad-blackoil ad-core ad-props mrst-gui

%% Set up grid
nx = 20;
ny = 20;
x_size = 1000;
y_size = 1000;
G = cartGrid([nx, ny], [x_size, y_size]);
G = computeGeometry(G);

%% Geometrical cells that are connected to the fractured domains
mrtr_cells = find((G.cells.centroids(:,1) >= 0.0 &...
                   G.cells.centroids(:,1) <= 1000));

%% Set up rock properties
base_rock = makeRock(G, 1*darcy, .01);
                
mrtr_matrix_rock = struct('perm',0.1*milli*darcy*ones(length(mrtr_cells),1),...
                                   'poro',.1*ones(length(mrtr_cells),1));

%% Set up fluid
base_fluid = initSimpleADIFluid('phases', 'WO', 'c', [1e-12,1e-12]/psia, 'n', [1,1],...
                                'mu',[1,1.5]*centi*poise);

sub_domains_fraction = [0.2; 0.35; 0.3; 0.1; 0.05];

mrtr_matrix_fluid = cell(length(sub_domains_fraction),1);
for i=1:length(sub_domains_fraction)
    mrtr_matrix_fluid{i} = initSimpleADIFluid('phases', 'WO', 'c', [1e-12,1e-12]/psia, 'n', [1,1],...
                                    'mu',[1,1.5]*centi*poise); 
end

%% Setting capillary pressure in the matrix of the dual porosity sub domains
for i=1:length(sub_domains_fraction)
    mrtr_matrix_fluid{i}.pcOW = @(sw)100*kilo*Pascal*sw.^-0.5;
end
%% Setting the porevolume multiplier to account for sub domain volume fractions
for i=1:length(sub_domains_fraction)
    mrtr_matrix_fluid{i}.pvMultR = @(p) p*0 + sub_domains_fraction(i);
end

%% Wells
W = [];
W = verticalWell(W, G, base_rock,  1,   1, (1:1),     ...
               'Type', 'rate', 'Val', 1e-4, ...
               'Name', 'P1', 'compi' ,[1,0]);     
W = verticalWell(W, G, base_rock,  nx,   ny, (1:1),     ...
               'Type', 'bhp', 'Val', 500*psia, ...
               'Name', 'P2', 'compi' ,[1,0]);    
           
%% Set up model
gravity reset off;
model = TwoPhaseOilWaterModel(G, base_rock, base_fluid);
model = model.validateModel();
                                                     
%% Fractured domain manager
manager = FracturedDomainManager();

% Adding matrix domains to the fracture cells with different constant transfer rates
transfer_model = cell(length(sub_domains_fraction),1);
transfer_model{1} = SaturationDifferenceTransferFunction(1e-8, 'maximum_saturation',1.0);
transfer_model{2} = SaturationDifferenceTransferFunction(1e-9, 'maximum_saturation',1.0);
transfer_model{3} = SaturationDifferenceTransferFunction(1e-10, 'maximum_saturation',1.0);
transfer_model{4} = SaturationDifferenceTransferFunction(1e-11, 'maximum_saturation',1.0);
transfer_model{5} = SaturationDifferenceTransferFunction(1e-12, 'maximum_saturation',1.0);

for i=1:length(sub_domains_fraction)
    model = manager.addFracturedDomain(model,...
                                         'multi_continuum',...
                                         mrtr_cells,...
                                         mrtr_matrix_rock,...
                                         mrtr_matrix_fluid{i},...
                                         'transfer_model',...
                                         transfer_model{i});     
end
                                 
%% Pore volume of the whole system (fracture cells plus matrix cells)

pv = sum(G.cells.volumes.*base_rock.poro);
for i=1:length(sub_domains_fraction)
    pv = pv + sub_domains_fraction(i)*sum(G.cells.volumes.*mrtr_matrix_rock.poro);
end

%% Initializing state 
eps = 1e-6;
state = initResSol(model.G, 1000*psia, [eps,1-eps]);
state.wellSol = initWellSolAD(W, model, state);
   
%% Solver
solver = NonLinearSolver();

%% Time loop
dt = 0.5*year;
tmax = 50*year;
t = 0;
time = (0);
oil_recovered=(0);
wcut=(0);

while t<=tmax
    
    disp(['Time = ',num2str(t/day), ' days'])
    state = solver.solveTimestep(state, dt, model, 'W', W);
    time = [time time(end) + dt];    
    oil_recovered = [oil_recovered oil_recovered(end) + abs(state.wellSol(2).qOs)*dt];
    wcut = [wcut state.wellSol(2).wcut];
    t = t+dt;
    
end

%% Plotting results
fig1 = figure('Position',[10,100,1400,700]);
fig1.Color = 'w';
fig2 = figure('Position',[500,100,400,400]);
fig2.Color = 'w';
fig3 = figure('Position',[900,100,400,400]);
fig3.Color = 'w';

% get the virtual cells from each domain for plotting purposes
vcells_m1 = model.G.FracturedDomains.domains{1}.virtual_cells;
vcells_m2 = model.G.FracturedDomains.domains{2}.virtual_cells;
vcells_m3 = model.G.FracturedDomains.domains{3}.virtual_cells;
vcells_m4 = model.G.FracturedDomains.domains{4}.virtual_cells;
vcells_m5 = model.G.FracturedDomains.domains{5}.virtual_cells;

% plot_cells are always geometrical cells
plot_cells = 1:nx*ny;

figure(fig1)
subplot(2,3,1)
colormap(flipud(jet));
p = plotCellData(G,state.s(plot_cells,1,1), plot_cells);
p.EdgeAlpha = 0.3;
colorbar;
caxis([0,1]);
set(gca,'FontSize',16);
xlabel('x')
ylabel('y')
hold on;
title('S_w (fracture cells)');

subplot(2,3,2)
colormap(flipud(jet));
p = plotCellData(G,state.s(vcells_m1,1,1), plot_cells);
p.EdgeAlpha = 0.3;
colorbar;
caxis([0,1]);
set(gca,'FontSize',16);
xlabel('x')
ylabel('y')
hold on;
title('S_w (matrix domain 1)');

subplot(2,3,3)
colormap(flipud(jet));
p = plotCellData(G,state.s(vcells_m2,1,1), plot_cells);
p.EdgeAlpha = 0.3;
colorbar;
caxis([0,1]);
set(gca,'FontSize',16);
xlabel('x')
ylabel('y')
hold on;
title('S_w (matrix domain 2)');

subplot(2,3,4)
colormap(flipud(jet));
p = plotCellData(G,state.s(vcells_m3,1,1), plot_cells);
p.EdgeAlpha = 0.3;
colorbar;
caxis([0,1]);
set(gca,'FontSize',16);
xlabel('x')
ylabel('y')
hold on;
title('S_w (matrix domain 3)');

subplot(2,3,5)
colormap(flipud(jet));
p = plotCellData(G,state.s(vcells_m4,1,1), plot_cells);
p.EdgeAlpha = 0.3;
colorbar;
caxis([0,1]);
set(gca,'FontSize',16);
xlabel('x')
ylabel('y')
hold on;
title('S_w (matrix domain 4)');

subplot(2,3,6)
colormap(flipud(jet));
p = plotCellData(G,state.s(vcells_m5,1,1), plot_cells);
p.EdgeAlpha = 0.3;
colorbar;
caxis([0,1]);
set(gca,'FontSize',16);
xlabel('x')
ylabel('y')
hold on;
title('S_w (matrix domain 5)');

figure(fig2)
plot(time/year,oil_recovered/pv, 'LineWidth', 2)
xlabel('Time [year]')
ylabel('Oil Recovery Factor [/]')
set(gca,'FontSize',16);

figure(fig3)
plot(time/year,wcut, 'LineWidth', 2)
xlabel('Time [year]')
ylabel('Water Cut [/]')
set(gca,'FontSize',16);

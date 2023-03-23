%% minc
clear;
clc;
close all;

mrstModule add ad-blackoil ad-core ad-props mrst-gui

%% Set up grid
nx = 3;
ny = 3;
x_size = 100;
y_size = 100;
G = cartGrid([nx, ny], [x_size, y_size]);
G = computeGeometry(G);

%% minc cells 
minc_cells = find((G.cells.centroids(:,1)>=0.0 &...
                            G.cells.centroids(:,1)<=1000));

%% Set up rock properties
base_rock = makeRock(G, 1000*milli*darcy, 100);
                
minc_matrix_rock = struct('perm',0.1*milli*darcy*ones(length(minc_cells),1),...
                                   'poro',.1*ones(length(minc_cells),1));
%% define the shell properties

shell_volume_fraction = [1/5, 1/5, 1/5, 1/5, 1/5];
% shell_volume_fraction = [1./27., 7./27., 19./27.];
% shell_volume_fraction = [1/2, 1/2];
% shell_volume_fraction = [0.8750, 0.0156];
% shell_volume_fraction = [1];
matrix_block_dimension = [1, 1, 1];

%% Set up fluid
base_fluid = initSimpleADIFluid('phases', 'WO', 'c', [1e-12,1e-12]/psia, 'n', [1,1],...
                                'mu',[1,1.5]*centi*poise);

minc_matrix_fluid = cell(length(shell_volume_fraction),1);
for i=1:length(shell_volume_fraction)
    minc_matrix_fluid{i} = initSimpleADIFluid('phases', 'WO', 'c', [1e-12,1e-12]/psia, 'n', [1,1],...
                                    'mu',[1,1.5]*centi*poise); 
end

%% Setting capillary pressure in the matrix of the dual porosity sub domains
for i=1:length(shell_volume_fraction)
    minc_matrix_fluid{i}.pcOW = @(sw)100*kilo*Pascal*sw.^-0.5;
end
%% Setting the porevolume multiplier to account for sub domain volume fractions
for i=1:length(shell_volume_fraction)
    minc_matrix_fluid{i}.pvMultR = @(p) p*0 + shell_volume_fraction(i);
end
%% wells
W = [];
           
%% Set up model
model = TwoPhaseOilWaterModel(G, base_rock, base_fluid);

%% Calling validate model
model = model.validateModel();

%% Fractured domain manager
manager = FracturedDomainManager();

transfer_model = cell(length(shell_volume_fraction),1);
for i=1:length(shell_volume_fraction)
transfer_model{i} = KazemiMultiphaseTransferFunction(MincShapeFactor(...
                    matrix_block_dimension,...
                    'volume_fractions', shell_volume_fraction,...
                    'domains', [i-1,i]));
end

% one can also define each flux between shells with constant rate
% transfer_model{1} = KazemiMultiphaseTransferFunction(...
%     ConstantShapeFactor(matrix_block_dimension, 72));
% transfer_model{2} = KazemiMultiphaseTransferFunction(...
%     ConstantShapeFactor(matrix_block_dimension, 16));
% transfer_model{3} = KazemiMultiphaseTransferFunction(...
%     ConstantShapeFactor(matrix_block_dimension, 8/3.));
% transfer_model{4} = KazemiMultiphaseTransferFunction(...
%     ConstantShapeFactor(matrix_block_dimension, 8/3.));
% transfer_model{5} = KazemiMultiphaseTransferFunction(...
%     ConstantShapeFactor(matrix_block_dimension, 8/3.));

minc_shell_cells = cell(length(shell_volume_fraction)+1,1);
minc_shell_cells{1} = minc_cells;
for i=1:length(shell_volume_fraction)
    model = manager.addFracturedDomain(model,...
                                         'multi_continuum',...
                                         minc_cells,...
                                         minc_matrix_rock,...
                                         minc_matrix_fluid{i},...
                                         'transfer_model',...
                                         transfer_model{i},...
                                         'connection_cell_list',...
                                         minc_shell_cells{i});     
    minc_shell_cells{i+1} = model.G.FracturedDomains.domains{i}.virtual_cells;
end

%% TO DO: solve problem
model.outputFluxes = 0;

%% Initializing state 
eps = 1e-6;
state = initResSol(model.G, 500*psia, [eps,1-eps]);
state.s(minc_cells,1)=1-eps;
state.s(minc_cells,2)=eps;
state.wellSol = initWellSolAD(W, model, state);
   
%% Solver
solver = NonLinearSolver();

%% Figure
% fig1 = figure('Position',[100,100,1200,800]);
% fig1.Color = 'w';
% colormap('jet');
fig2 = figure('Position',[100,100,600,600]);
fig2.Color = 'w';

[ dist, area] = shell_properties( shell_volume_fraction, matrix_block_dimension )
x_vct = [0.0, dist(1)];
for i=1:length(dist)-1
    x_vct = [x_vct x_vct(end) + dist(i)+dist(i+1)];
end

%% Time loop
dt = 1*day;
tmax = 20*dt;
t = 0;
while t<=tmax
    
    disp(['Time = ',num2str(t/day), ' days'])
    state = solver.solveTimestep(state, dt, model, 'W', W);
    state.s(minc_cells,1)=1-eps;
    state.s(minc_cells,2)=eps;

    
%     figure(fig1)
%     subplot(4,2,1);
%     p = plotCellData(G,state.s(minc_shell_cells{1},1,1));
%     p.EdgeAlpha = 0;
%     colorbar;
%     caxis([0,1]);
%     set(gca,'FontSize',16);
%     xlabel('x')
%     ylabel('y')
%     hold on;
%     
%     figure(fig1)
%     subplot(4,2,2);
%     p = plotCellData(G,state.s(minc_shell_cells{2},1),minc_cells);
%     p.EdgeAlpha = 0;
%     colorbar;
%     caxis([0,1]);
%     xlim([0,x_size]);
%     ylim([0,y_size]);
%     set(gca,'FontSize',16);
%     xlabel('x')
%     ylabel('y')
%     drawnow;
%     
%     figure(fig1)
%     subplot(4,2,3);
%     p = plotCellData(G,state.s(minc_shell_cells{3},1),minc_cells);
%     p.EdgeAlpha = 0;
%     colorbar;
%     caxis([0,1]);
%     xlim([0,x_size]);
%     ylim([0,y_size]);
%     set(gca,'FontSize',16);
%     xlabel('x')
%     ylabel('y')
%     drawnow;
% 
%     figure(fig1)
%     subplot(4,2,4);
%     p = plotCellData(G,state.s(minc_shell_cells{4},1),minc_cells);
%     p.EdgeAlpha = 0;
%     colorbar;
%     caxis([0,1]);
%     xlim([0,x_size]);
%     ylim([0,y_size]);
%     set(gca,'FontSize',16);
%     xlabel('x')
%     ylabel('y')
%     drawnow;
% 
%     figure(fig1)
%     subplot(4,2,5);
%     p = plotCellData(G,state.pressure(minc_shell_cells{1})/psia);
%     p.EdgeAlpha = 0;
%     colorbar;
% %     caxis([0,1]);
%     set(gca,'FontSize',16);
%     xlabel('x')
%     ylabel('y')
%     hold on;
%     
%     figure(fig1)
%     subplot(4,2,6);
%     p = plotCellData(G,state.pressure(minc_shell_cells{2})/psia,minc_cells);
%     p.EdgeAlpha = 0;
%     colorbar;
% %     caxis([0,1]);
%     xlim([0,x_size]);
%     ylim([0,y_size]);
%     set(gca,'FontSize',16);
%     xlabel('x')
%     ylabel('y')
%     drawnow;
%     
%     figure(fig1)
%     subplot(4,2,7);
%     p = plotCellData(G,state.pressure(minc_shell_cells{3})/psia,minc_cells);
%     p.EdgeAlpha = 0;
%     colorbar;
% %     caxis([0,1]);
%     xlim([0,x_size]);
%     ylim([0,y_size]);
%     set(gca,'FontSize',16);
%     xlabel('x')
%     ylabel('y')
%     drawnow;
% 
%     figure(fig1)
%     subplot(4,2,8);
%     p = plotCellData(G,state.pressure(minc_shell_cells{4})/psia,minc_cells);
%     p.EdgeAlpha = 0;
%     colorbar;
% %     caxis([0,1]);
%     xlim([0,x_size]);
%     ylim([0,y_size]);
%     set(gca,'FontSize',16);
%     xlabel('x')
%     ylabel('y')
%     drawnow;

    sw = [];
    for i=1:length(minc_shell_cells)
        sw = [sw state.s(minc_shell_cells{i}(5),1)];
    end
    pressure = [];
    for i=1:length(minc_shell_cells)
        pressure = [pressure state.pressure(minc_shell_cells{i}(5))/psia];
    end
    
    figure(fig2)
    grid on;
    subplot(2,1,1);
    hold on;
    p = plot(x_vct,sw,'DisplayName',num2str(t/day),'Marker','o');
%     p = plot(state.s(minc_shell_cells{1},1),'DisplayName','fracture');
%     hold on;
%     p = plot(state.s(minc_shell_cells{2},1),'DisplayName','shell1');
%     p = plot(state.s(minc_shell_cells{3},1),'DisplayName','shell2');
%     p = plot(state.s(minc_shell_cells{4},1),'DisplayName','shell3');
    set(gca,'FontSize',16);
    xlabel('x')
    ylabel('saturation')
%     hold off;
    drawnow;
    
    figure(fig2)
    grid on;
    subplot(2,1,2);
    hold on;
    p = plot(x_vct,pressure,'DisplayName',num2str(t/day),'Marker','o');
%     p = plot(state.pressure(minc_shell_cells{1})/psia,'DisplayName','fracture');
%     hold on;
%     p = plot(state.pressure(minc_shell_cells{2})/psia,'DisplayName','shell1');
%     p = plot(state.pressure(minc_shell_cells{3})/psia,'DisplayName','shell2');
%     p = plot(state.pressure(minc_shell_cells{4})/psia,'DisplayName','shell3');
    set(gca,'FontSize',16);
    xlabel('x')
    ylabel('pressure')
%     hold off;
    drawnow;

    t = t+dt;
    
end


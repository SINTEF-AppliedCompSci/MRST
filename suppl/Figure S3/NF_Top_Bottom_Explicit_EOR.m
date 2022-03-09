%% Loading Modules
clear; 
clc;
close all;
Globals;
case2run ='ProdBot_InjTop';%'ProdBot_InjTop'; %'ProdBot';
opt = struct('nkr',        1, ...
             'shouldPlot', 1 ); %change to 0 if running on HPC
 %% Load necessary modules, etc 
mrstModule add hfm;             % hybrid fracture module       
mrstModule add ad-props ad-core 
mrstModule add mrst-gui;        % plotting routines
mrstModule add compositional; %mrstModule add ad-blackoil;
mrstModule add upr;
%% Basic Parameters
tol=1e-5;    
physdim = [200 100 80]*meter; %sets grid's physical dimensions (x,y,z)
wf = 10*milli*meter; %SlotFrac aperture
fracloc = [10 70]*meter; %location of SlotFrac
layer_z = 20; %No. of layers in z-direction, excluding Slotfrac layers
layer_frac = fix(layer_z/(numel(fracloc)*2)); %calculate no. of layers above or under frac plan   
%% Create a Box Reservoir Grid Using Geometric Grid (OMO):
% Z gridding for SD Fractures
% space on top of top frac
endGeom = 10;
numPts = 20; %20
dist2grid = fracloc(1);
ptZ = geomspace(wf/2, dist2grid, numPts, wf);
refZ = fracloc(1);
facesZcoords_aboveTop = refZ - fliplr(ptZ);
facesZcoords_aboveTop(1)=0;
% space below bottom frac
refZ = fracloc(2);
facesZcoords_belowBottom = refZ + ptZ;
% space above bottom frac
numPts = 15; %15
dist2grid = endGeom;%0.5*sum(fracloc)-fracloc(1); 
ptZ = geomspace(wf/2, dist2grid, numPts, wf);
facesZcoords_aboveBottom = refZ - fliplr(ptZ);
% space below top frac
refZ = fracloc(1);
facesZcoords_belowTop = refZ + ptZ;
% space in-between two SD fractures starting from endGeom point 1 and 2
last_log_interval = facesZcoords_belowTop(end)-facesZcoords_belowTop(end-1);
linear_points = ceil(((fracloc(2)-endGeom)-(fracloc(1)+endGeom))/last_log_interval);
facesZcoords_between = linspace(fracloc(1)+endGeom,fracloc(2)-endGeom,linear_points);
z = [facesZcoords_aboveTop, facesZcoords_belowTop(1:end-1),...
                facesZcoords_between,facesZcoords_aboveBottom(2:end), facesZcoords_belowBottom];

% X-gridding for NF
[G_NF,frac_ids,well_ids]=explicitStencilSuppl(physdim, 'aperture',0.1*ft,...
                             'numStages',1,'fracSpacing', 200, 'fracHalfLength'  , 39.99,      ...
                             'nxRefine',64,'tipNY',2,'ny',8, 'nz',8, 'gridType','geomspacing'); 

G_matrix = tensorGrid(G_NF.XGrid, 0:5:physdim(2),z);
% G_matrix = tensorGrid(0:5:physdim(1), 0:5:physdim(2),z);
G_matrix = computeGeometry(G_matrix);    
G_matrix.rock=makeShaleRock(G_matrix,10*micro*darcy,0.07);
    
frac_z=[];
for ii = 1:length(fracloc) %this nested loop locate fracture z-layer index from fracture location.
    n_layer = 0;
    for jj = 1:length(z)-1
        if ~((fracloc(ii) > z(jj)) && (fracloc(ii) < z(jj+1)))
            n_layer = n_layer+1;
        else
            frac_z = [frac_z,n_layer+1];
            break;
        end
    end
end
% %Model SD fracture explicitly using a semicircle (only applied for an x*y*z
% %reservoir where x = 2y)
% fraccells=[];
% element_size = physdim ./ G_matrix.cartDims;
% left_end = physdim(1)/2;
% for i = 1:numel(frac_z)
%     for j = 1:G_matrix.cartDims(1)/2
%         temp = 0; jj =0 ;
%         Frac_I = j;
%         Frac_I_flip = G_matrix.cartDims(1)- j;
%         idx_start = left_end - j*element_size(1);
%         while temp <= sqrt(left_end^2-idx_start^2) %+element_size(2)
%             jj = jj +1;
%             temp = jj*element_size(2);
%         end
%         Frac_J = G_matrix.cartDims(2):-1:G_matrix.cartDims(2)-jj+2;
%         Frac_J_flip = G_matrix.cartDims(2):-1:jj+2;
%         Frac_K = frac_z(i);
%         [II,JJ,KK]=meshgrid(Frac_I,Frac_J,Frac_K);
%         [II_flip,JJ_flip,KK_flip]=meshgrid(Frac_I_flip,Frac_J,Frac_K);
% %         [II_flip,JJ_flip,KK_flip]=meshgrid(Frac_I_flip,Frac_J_flip,Frac_K);
%         FracCellIds = sub2ind(G_matrix.cartDims, II,JJ,KK);
%         FracCellIds_flip = sub2ind(G_matrix.cartDims, II_flip,JJ_flip,KK_flip);
%         fraccells=[fraccells; FracCellIds(:);FracCellIds_flip(:)];
%     end
% end
[fracIndx_X,fracIndx_Y,fracIndx_Z] = meshgrid(1:G_matrix.cartDims(1), 1:G_matrix.cartDims(2), frac_z);
fraccells = sub2ind(G_matrix.cartDims, reshape(fracIndx_X,numel(fracIndx_X),1),reshape(fracIndx_Y,...
            numel(fracIndx_X),1), reshape(fracIndx_Z,numel(fracIndx_X),1));
G_matrix.rock.poro(fraccells) = 0.33;
G_matrix.rock.perm(fraccells) = 10*darcy; 

%Assign NF parameters
[NFIndx_X,NFIndx_Y,NFIndx_Z] = meshgrid(G_NF.Frac_I, 4:G_matrix.cartDims(2)-4, 4:G_matrix.cartDims(3)-3);
NFcells = sub2ind(G_matrix.cartDims, reshape(NFIndx_X,numel(NFIndx_X),1),reshape(NFIndx_Y,...
            numel(NFIndx_X),1), reshape(NFIndx_Z,numel(NFIndx_X),1));
G_matrix.rock.poro(NFcells) = 0.5;
G_matrix.rock.perm(NFcells) = 10*milli*darcy; 

G=G_matrix;

if (opt.shouldPlot)
    figure
    %Plot matrix
    show = true([G.cells.num, 1]);
    show(fraccells(:)) = false;% Hide well cell
    plotCellData (G , convertTo(G.rock.perm,milli*darcy),show, ...
        'EdgeColor', 'k','facealpha',0.0);
    colorbar ('horiz'); view(40,30); axis equal tight;
    hold on;
    %Plot frac plane
    show = false([G.cells.num, 1]);
    show(fraccells(:)) = true;% Hide well cell
    show(NFcells(:)) = true;
    plotCellData (G , convertTo(G.rock.perm,milli*darcy),show, ...
        'EdgeColor', 'k');

    hold off;
end
%% Define three-phase compressible flow model
useNatural = true;
casename = 'oil_1';
pwf = 2500*psia; 
rate = 0.003277; %10,000 scf/day = 0.003277 m^3/s

[fluid, info] = getShaleCompFluidCase(casename);
eosname = 'prcorr';  %'srk','rk','prcorr'
G1cell = cartGrid([1 1],[1 1]);
G1cell = computeGeometry(G1cell);
EOSModel = EquationOfStateModel(G1cell, fluid, eosname);

%Surface Conditions
p_sc = 101325; %atmospheric pressure
T_sc = 288.706;% 60 Farenheit
[L, x, y, Z_L, Z_V, rhoO_S, rhoG_S] = standaloneFlash(p_sc, T_sc, info.initial, EOSModel);
flowfluid = initSimpleADIFluid('phases', 'WOG', 'n', [opt.nkr, opt.nkr, opt.nkr], 'rho', [1000, rhoO_S, rhoG_S]);    % flowfluid.KGangiFn = @(p) power((1-power(((Pc - alpha.*p)./Pmax),m)),3);
gravity reset on

arg = {G, G.rock, flowfluid, fluid, 'water', true};

diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);
mex_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'useMex', true);
sparse_backend = SparseAutoDiffBackend();

if useNatural
    constructor = @NatVarsShaleModel;
else
    constructor = @GenericOverallCompositionModel;
end

modelSparseAD = constructor(arg{:}, 'AutoDiffBackend', sparse_backend);
modelDiagonalAD = constructor(arg{:}, 'AutoDiffBackend', diagonal_backend);
modelMexDiagonalAD = constructor(arg{:}, 'AutoDiffBackend', mex_backend);

% modelSparseAD.operators = TPFAoperators;
% modelDiagonalAD.operators = TPFAoperators;
% modelMexDiagonalAD.operators = TPFAoperators;
%% Set up initial state
totTime = 15*year;
nSteps =15;
ncomp = fluid.getNumberOfComponents();
s0 = [0.23, 0.70, 0.07];  
state = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, modelSparseAD.EOSModel);
%% Pick linear solver
% The AMGCL library is one possible solver option for MRST. It is fairly
% easy to write interfaces to other solvers using MEX files and/or the
% LinearSolverAD class. We call the routine for automatically selecting a
% reasonable linear solver for the specific model.
linsolve = selectLinearSolverAD(modelDiagonalAD,'useAMGCL',true,'useCPR',true);
disp(linsolve)
nls = NonLinearSolver('LinearSolver', linsolve);
%% Schedule
wellRadius = 0.01;
W = [];
switch case2run
case 'ProdTop' 
    % Producer
    W = verticalWell(W, G.Matrix, G.Matrix.rock, 1, 1, 10, ...
        'comp_i', [0.23, 0.76, 0.01],'Name', 'Prod_Top', 'Val', pwf, 'sign', -1, 'Type', 'bhp','Radius', wellRadius);
    W(1).components = info.initial;
case 'ProdTop_InjBot'
    % Producer
    W = verticalWell(W, G.Matrix, G.Matrix.rock, 2, 2, 10, ...
        'comp_i', [0.23, 0.76, 0.01],'Name', 'Prod_Top', 'Val', pwf, 'sign', -1, 'Type', 'bhp','Radius', wellRadius);
    % Injector
    W = verticalWell(W, G.Matrix, G.Matrix.rock, 1, 1, 70, ...
        'comp_i', [0 0 1],'Name', 'Inj_Bot', 'Val', pinj, 'sign', 1, 'Type', 'bhp','Radius', wellRadius);
    W(1).components = info.initial;
    W(2).components = info.injection;
case 'ProdBot' 
    % Producer
    W = verticalWell(W, G.Matrix, G.Matrix.rock, 5, 3, frac_z(2), ...
        'comp_i', [0.23, 0.76, 0.01],'Name', 'Prod_Bot', 'Val', pwf, 'sign', -1, 'Type', 'bhp','Radius', wellRadius);
    W(1).components = info.initial;
case 'ProdBot_InjTop'
    % Producer
    W = verticalWell(W, G, G.rock, 20, 10, frac_z(2), ...
        'comp_i', [0.23, 0.76, 0.01],'Name', 'Prod_Bot', 'Val', pwf, 'sign', -1, 'Type', 'bhp','Radius', wellRadius); 
%     W = verticalWell(W, G.Matrix, G.Matrix.rock, 20, 10, frac_z(1), ...
%         'comp_i', [0 0 1],'Name', 'Inj_Top', 'Val', pinj, 'sign', 1, 'Type', 'bhp','Radius', wellRadius); %control by injection pressure
    % Injector
    W = verticalWell(W, G, G.rock, 20, 10, frac_z(1), ...
        'comp_i', [0 0 1],'Name', 'Inj_Top', 'Val', rate, 'sign', 1, 'Type', 'rate','Radius', wellRadius); %control by injection rate
    W(1).components = info.initial;
    W(2).components = info.injection;
otherwise
    warning('Case Does Not Exist. Running case with Prod Only at Bottom')
    % Producer
    W = verticalWell(W, G.Matrix, G.Matrix.rock, 1, 1, 80, ...
        'comp_i', [0.23, 0.76, 0.01],'Name', 'Prod_Bot', 'Val', pwf, 'sign', -1, 'Type', 'bhp','Radius', wellRadius);
    W(1).components = info.initial;
end
dt = rampupTimesteps(totTime, 30*day, nSteps); %20*day
schedule = simpleSchedule(dt, 'W', W);  
%% Simulate problem
% frac_intensity = sum(fracArea)/prod(physdim,'all')
% [ws, states, reports] = simulateScheduleAD(state, modelMexDiagonalAD, schedule, 'nonlinearsolver', nls, 'Verbose', true);
% [ws, states, reports] = simulateScheduleAD(state, modelDiagonalAD, schedule, 'nonlinearsolver', nls, 'Verbose', true);
[ws, states, reports] = simulateScheduleAD(state, modelMexDiagonalAD, schedule, 'nonlinearsolver', nls, 'Verbose', true);
% [ws, states, reports] = simulateScheduleAD(state, modelDiagonalAD, schedule, 'Verbose', true); %direct solver
%% plotting
if (opt.shouldPlot)
    figure, 
    pargs = {'EdgeColor','k'};
    plotToolbar(G, states,pargs{:})
    view(40,30);
    axis tight equal;
    plotWellSols(ws,cumsum(schedule.step.val))
end
% Calculate RF
tinSecs = cumsum(schedule.step.val);
tinDays = tinSecs./86400;
Boi = ws{1,1}(1).qOr/ws{1,1}(1).qOs; %calculate the initial oil formation volume factor from production rates
STOIIP = 6.289811*sum((G_matrix.cells.volumes .* G_matrix.rock.poro))*(1-s0(1))/Boi; %in STB

qO = [];
for ii=1:size(ws)
    qO = [qO,-543439.65*ws{ii,1}(1).qOs]; %covnert from m^3/s to stb/d
end
Np = trapz(tinDays,qO);
RF = Np/STOIIP;
%% Get cumulative reservoir fluid withdrawn
x = zeros(length(ws),1);
for i = 1:length(ws)
   x(i)= -543439.7*ws{i}(1).qTr;
end
QTr = trapz(tinDays,x);
%% Save Output Variables (Used in HPC).
% if ~opt.shouldPlot
%     fpath =  '/scratch/ahass16/';
%     fullFinalOut = [fpath, 'NF_Top_Bottom_Explicit_EOR.mat'];
%     save(fullFinalOut,'-v7.3');
% end
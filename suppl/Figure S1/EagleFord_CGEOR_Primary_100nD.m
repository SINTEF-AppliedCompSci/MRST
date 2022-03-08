%% Loading Modules
clear; 
clc;
close all;
Globals;
opt = struct('nkr',        1, ...
             'shouldPlot', 0 );
%% Load necessary modules, etc 
mrstModule add hfm;             % hybrid fracture module  
mrstModule add ad-props ad-core 
mrstModule add mrst-gui;        % plotting routines
mrstModule add compositional; %mrstModule add ad-blackoil;
mrstModule add upr;
%% Basic Parameters
tol=1e-5;    
physdim = [200 100 80]*meter; %sets grid's physical dimensions (x,y,z)
fracHalfLength=30;
fracHeight=40;
fractureSpacing=37.5;
fracPerStage=2;
NumStages=3;
clusterSpacing=15;
fracAperture=0.0030;
fracPoro=0.5;
fracPerm=1.0*darcy;
MatrixPoro=0.07;
MatrixPerm=100*nano*darcy;
wellLength=140;
wellRadius=0.1000;
heelCooord=[30, physdim(2)/2, physdim(3)/2];
EndCoord=[heelCooord(1)+wellLength,physdim(2)/2, physdim(3)/2];

heelCooord_inj=[30+1e-5, physdim(2)/2, physdim(3)/2];
EndCoord_inj=[heelCooord(1)+wellLength,physdim(2)/2, physdim(3)/2];
%% Create Fracture System
[fracplanes,frac_centroid_s] = createMultiStageHFsCGEOR('numStages',NumStages,'fracSpacing', fractureSpacing,...
      'numFracsPerStage', fracPerStage,'fracHalfLength', fracHalfLength,'fracHeight',fracHeight, ...
      'clusterSpacing', clusterSpacing,'heelCoord',[40, physdim(2)/2, physdim(3)/2]);
 
for i=1:numel(fracplanes)
     fracplanes(i).aperture = fracAperture; 
     fracplanes(i).poro = fracPoro;
     fracplanes(i).perm = fracPerm;
end
G_matrix = meshHFsystem(physdim,frac_centroid_s,'numStages',NumStages,'numFracsPerStage',fracPerStage,...
    'ny',21,'nz',17,'nxRefine_small',4,'nxRefine',10,'fracSpacing', fractureSpacing,'aperture',fracAperture);
G_matrix = computeGeometry(G_matrix);
G_matrix.rock=makeRock(G_matrix,MatrixPerm,MatrixPoro);

if (opt.shouldPlot)
    plotfracongrid(G_matrix,fracplanes,'label',false); % visualize to check before pre-process
    view(40,30); 
end
%% Create Wells
wells = struct;
% wells(1).points=[heelCooord_inj; EndCoord_inj];
wells(1).points=[heelCooord; EndCoord];
wells(1).radius=wellRadius;
wells(2).points=[heelCooord; EndCoord];
wells(2).radius=wellRadius;
%% EDFM PreProcessing
G=G_matrix;
[G,fracplanes]=EDFMshalegrid(G,fracplanes,...
        'Tolerance',tol,'plotgrid',false,'fracturelist',1:numel(fracplanes));
%-Frac-Matrix NNCs
G = fracturematrixShaleNNC3D(G,tol);
%-Frac-Frac NNCs
[G,fracplanes]=fracturefractureShaleNNCs3D(G,fracplanes,tol);
%-Well-Fracs NNCs
[G,wells] = wellfractureNNCs3D(G,fracplanes,wells,tol);

if (opt.shouldPlot)
    figure,
    plotfracSystem_NF(G,fracplanes,numel(fracplanes),wells,'label',false)
    view(40,30);
end
TPFAoperators = setupShaleEDFMOpsTPFA(G, G.rock, tol);
%% Define three-phase compressible flow model
useNatural = true;
casename = 'eagleford';
pwf = 2500*psia;
rate  = 0.03277; %1e5
[fluid, info] = getShaleCompFluidCase(casename);
info.pressure = 7000*psia;
% info.injection = [0 0 1 0 0 0 0]; %C1 injection
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
arg = {G, [], flowfluid, fluid, 'water', true};

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

modelSparseAD.operators = TPFAoperators;
modelDiagonalAD.operators = TPFAoperators;
modelMexDiagonalAD.operators = TPFAoperators;
%% Set up initial state
ncomp = fluid.getNumberOfComponents();
s0 = [0.23, 0.70, 0.07]; 
state = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, modelSparseAD.EOSModel);
%% Pick linear solver
linsolve = selectLinearSolverAD(modelDiagonalAD,'useAMGCL',true,'useCPR',true);
disp(linsolve)
nls = NonLinearSolver('LinearSolver', linsolve);
%% HnP Well Definition
W = [];
W = addWellEDFMshale(W, G, G.Matrix.rock, wells(1).XFracCellIDs, ...
    'comp_i', [0, 0, 1],'Name', 'Injector', 'Val',...
    rate, 'sign', 1, 'Type', 'rate','Radius', wells(1).radius,'Dir','x'); %injector

W = addWellEDFMshale(W, G, G.Matrix.rock, wells(2).XFracCellIDs, ...
    'comp_i', [0.17, 0.74, 0.09],'Name', 'Producer', 'Val',...
    pwf, 'sign', -1, 'Type', 'bhp','Radius', wells(2).radius,'Dir','x'); %Producer
W(1).components = info.injection;
W(2).components = info.initial;
%% HnP Schedule
time_inj = [60*day 10*day 15]; %
time_soak = [14*day 3*day 5];
time_prod = [180*day 40*day 5]; %DOE/NTEL White Paper IOR Report
primaryTime = [3*year 30*day 15]; %totTime, 30*day, nSteps
HnPTime = 8*year;
schedule = Primary_HnP_schedule(primaryTime,time_inj,time_soak,time_prod,HnPTime, W);
tinSecs = cumsum(schedule.step.val);
tinDays = tinSecs./86400;
%% Simulate HnP schedule
[ws, states, reports] = simulateScheduleAD(state, modelSparseAD, schedule, 'nonlinearsolver', nls, 'Verbose', true);
%% plotting
if (opt.shouldPlot)
    figure, 
    plotToolbar(G, states)
    view(40,30);
    axis tight equal;
    plotWellSols(ws,cumsum(schedule.step.val))
%     plotWell(G,W); 
end
%% Calculate RF
tinSecs = cumsum(schedule.step.val);
tinDays = tinSecs./86400;
% Boi = ws{1,1}(1).qOr/ws{1,1}(1).qOs; %calculate the initial oil formation volume factor from production rates
% STOIIP = 6.289811*sum((G_matrix.cells.volumes .* G_matrix.rock.poro))*(1-s0(1))/Boi; %in STB
% 
% qO = [];
% for ii=1:size(ws)
%     qO = [qO,-543439.65*ws{ii,1}(1).qOs]; %covnert from m^3/s to stb/d
% end
% Np = trapz(tinDays,qO);
% RF = Np/STOIIP;
% %% Get cumulative reservoir fluid withdrawn
% x = zeros(length(ws),1);
% for i = 1:length(ws)
%    x(i)= -543439.7*ws{i}(1).qTr;
% end
% QTr = trapz(tinDays,x);
%% Save Output Variables (Used in HPC).
if ~opt.shouldPlot
    fpath =  '/scratch/ahass16/';
    fullFinalOut = [fpath, 'EagleFord_CGEOR_Primary_100nD.mat'];
    save(fullFinalOut,'ws','G','schedule');
end
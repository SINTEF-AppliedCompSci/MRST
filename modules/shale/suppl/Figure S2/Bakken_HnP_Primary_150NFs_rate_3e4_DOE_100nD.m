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
MatrixPerm=100*nano*darcy; %10*micro*darcy; %100*nano*darcy
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
%% Microfractures parallel to bedding planes
% rng(123456789); 
% tol4domain = tol*1e3;
% 
% set1 = Field(DFN('dim',3,'n',50,'dir',100,'ddir',0,'minl',5,...
%            'mu',10,'maxl',15,'bbx',[tol4domain,tol4domain,tol4domain,physdim(1)-tol4domain,physdim(2)-tol4domain,physdim(3)-tol4domain],'dip',-45,'ddip',0,...
%            'shape','s'),'Poly'); 
% set2 = Field(DFN('dim',3,'n',50,'dir',45,'ddir',-1e9,'minl',5,...
%             'mu',10,'maxl',15,'bbx',[tol4domain,tol4domain,tol4domain,physdim(1)-tol4domain,physdim(2)-tol4domain,physdim(3)-tol4domain],'dip',45,'ddip',-1e9,...
%             'shape','s'),'Poly');  
% set3 = Field(DFN('dim',3,'n',50,'dir',-55,'ddir',-1e9,'minl',5,...
%             'mu',10,'maxl',15,'bbx',[tol4domain,tol4domain,tol4domain,physdim(1)-tol4domain,physdim(2)-tol4domain,physdim(3)-tol4domain],'dip',135,'ddip',-1e9,...
%             'shape','s'),'Poly'); %'l','q',4

rng(123456789);  %(simple NF)
tol4domain = tol*1e3;

set1 = Field(DFN('dim',3,'n',4,'dir',100,'ddir',0,'minl',5,...
           'mu',10,'maxl',15,'bbx',[tol4domain,tol4domain,tol4domain,physdim(1)-tol4domain,physdim(2)-tol4domain,physdim(3)-tol4domain],'dip',-45,'ddip',0,...
           'shape','s'),'Poly'); 
set2 = Field(DFN('dim',3,'n',4,'dir',45,'ddir',-1e9,'minl',5,...
            'mu',10,'maxl',15,'bbx',[tol4domain,tol4domain,tol4domain,physdim(1)-tol4domain,physdim(2)-tol4domain,physdim(3)-tol4domain],'dip',45,'ddip',-1e9,...
            'shape','s'),'Poly');  
set3 = Field(DFN('dim',3,'n',4,'dir',-55,'ddir',-1e9,'minl',5,...
            'mu',10,'maxl',15,'bbx',[tol4domain,tol4domain,tol4domain,physdim(1)-tol4domain,physdim(2)-tol4domain,physdim(3)-tol4domain],'dip',135,'ddip',-1e9,...
            'shape','s'),'Poly'); %'l','q',4
			
set2 = [set2;set3];
[set1_,nonPlanarSets1,fracArea1] = processStochFracs(set1);
[set2_,nonPlanarSets2,fracArea2] = processStochFracs(set2); 
fracArea = [fracArea1;fracArea2];
% frac_intensity = sum(fracArea)/prod(physdim,'all')
fprintf('%d of %d set1 fracs were OK while %d of %d set2 fracs were OK \n',...
    numel(set1_),numel(set1),numel(set2_),numel(set2));

fprintf('Number of nonplanar sets in sets 1 and 2 are : %d and %d respectively\n',...
    numel(nonPlanarSets1),numel(nonPlanarSets2));
if (opt.shouldPlot)
    figure;
    Draw('ply',set1_);
    Draw('ply',set2_);view(45,30)
end 
fracSet = [set1_ ;set2_]; 
numHFplanes = numel(fracplanes);
numNFplanes = numel(fracSet); 
totalNumFracPlanes = numHFplanes + numNFplanes;
startIdx = numHFplanes + 1;
for i=1:numNFplanes
    idxGlobal = numHFplanes + i;
    fracplanes(idxGlobal).points = fracSet{i}(1:end-1,:);
    fracplanes(idxGlobal).aperture =  0.1*ft;
    fracplanes(idxGlobal).poro=0.5;
    fracplanes(idxGlobal).perm=10*milli*darcy;
end
G_matrix.numHFplanes=numHFplanes;
G_matrix.numNFplanes=numNFplanes;

if (opt.shouldPlot)
    figure,
    plotfracongrid(G_matrix,fracplanes,'label',false); % visualize to check before pre-process
end
G=G_matrix;
%% EDFM PreProcessing
G.rock.shaleMechanisms.Diffusion = false;
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
casename = 'bakken_light_5comps';
pwf = 3000*psia;
rate  = 0.003277*3; %1e4

[fluid, info] = getShaleCompFluidCase(casename);
info.pressure = 6700*psia;
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
W = addWellShale(W, G, G.Matrix.rock, wells(1).XFracCellIDs, ...
    'comp_i', [0, 0, 1],'Name', 'Injector', 'Val',...
    rate, 'sign', 1, 'Type', 'rate','Radius', wells(1).radius,'Dir','x'); %injector

W = addWellShale(W, G, G.Matrix.rock, wells(2).XFracCellIDs, ...
    'comp_i', [0.17, 0.74, 0.09],'Name', 'Producer', 'Val',...
    pwf, 'sign', -1, 'Type', 'bhp','Radius', wells(2).radius,'Dir','x'); %Producer
W(1).components = info.injection;
W(2).components = info.initial;
%% HnP Schedule
time_inj = [60*day 10*day 12]; %
time_soak = [14*day 3*day 5];
time_prod = [180*day 40*day 5]; %DOE/NTEL White Paper IOR Report
primaryTime = [3*year 90*day 15]; %totTime, 30*day, nSteps
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
% 
% data1(isnan(data1))=0;
% data11(isnan(data11))=0;
% qos = cumtrapz(tinDays,data1)./1e3;
% qgs = cumtrapz(tinDays,data11)./1e6;
% tinyrs = tinDays./365.25;
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
    fullFinalOut = [fpath, 'Bakken_HnP_Primary_NFs_rate_3e4_DOE.mat'];
    save(fullFinalOut,'ws','G','schedule');
end
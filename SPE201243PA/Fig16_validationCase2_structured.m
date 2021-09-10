%% This will generate Figure 16 in SPE-201243-PA
%
clear 
clc
close all
opt = struct('nkr',        2, ...
    'shouldPlot', 1 );

%% Load necessary modules, etc
mrstModule add hfm;             % hybrid fracture module
mrstModule add ad-core;
mrstModule add ad-props ad-core
mrstModule add mrst-gui;        % plotting routines
mrstModule add compositional; %mrstModule add ad-blackoil;
mrstModule add upr;

%% Basic Parameters
tol=1e-3;
physdim = [1990 1990 150]*ft;

wellRadius = 0.25*ft;

matrixPoro = 0.16;
matrixPerm = 100*micro*darcy;

w_f = 0.01*ft;
k_f = 50*darcy;
x_f = 350*ft;


%% Create Box Reservoir Grid

% [G_matrix,~,~]=explicitStencil3Fracs(physdim, 'aperture',w_f,...
%       'numStages',1,'fracSpacing',physdim(1), 'fracHalfLength'  , x_f,      ...
%        'nxRefine',12,'nxRefineSmall',10,'tipNY',8,'ny',8, 'nz',6, 'gridType','cartesian');    
[G_matrix,frac_ids,well_ids]=explicitStencil(physdim, 'aperture',w_f,...
      'numStages',1,'fracSpacing', physdim(1), 'fracHalfLength'  , x_f,...
       'nxRefine',32,'tipNY',8,'ny',8, 'nz',6, 'gridType','cartesian'); 
G_matrix = computeGeometry(G_matrix);

if (opt.shouldPlot)
    figure,
    plotGrid(G_matrix)%, view(5,45)
end

G_matrix.rock = makeShaleRock(G_matrix, matrixPerm, matrixPoro);

%% Create Fracture System

 %-Hydraulic Fracture Creation
fracplanes = struct;

fracplanes(1).points = [930 645 0; 930 1345 0; 930 1345 150; 930 645 150]*ft;
fracplanes(2).points = [995 645 0; 995 1345 0; 995 1345 150; 995 645 150]*ft;
fracplanes(3).points = [1060 645 0; 1060 1345 0; 1060 1345 150; 1060 645 150]*ft;

barrierPerm = 1*nano*darcy;
%barrier 1
fracplanes(1).aperture = w_f;
fracplanes(1).poro= matrixPoro;
fracplanes(1).perm= barrierPerm;

%fracture
fracplanes(2).aperture = w_f;
fracplanes(2).poro= matrixPoro;
fracplanes(2).perm= k_f;

%barrier 2
fracplanes(3).aperture = w_f;
fracplanes(3).poro= matrixPoro;
fracplanes(3).perm= barrierPerm;

if (opt.shouldPlot)
    plotfracongrid(G_matrix,fracplanes); % visualize to check before pre-process
end

%% Create Wells

wells = struct;
%We use a short horizontal well to allow us find fracture/well intersection
halfToe2HeelDist = 5;
midpoint = physdim/2;
wells(1).points=[midpoint(1)-halfToe2HeelDist,midpoint(2),midpoint(3); midpoint(1)+halfToe2HeelDist,midpoint(2),midpoint(3)]*meter;
wells(1).radius=wellRadius;

%% visualize to check before pre-process
G=G_matrix;

%% EDFM PreProcessing
[G,fracplanes]=EDFMshalegrid(G,fracplanes,...
        'Tolerance',tol,'plotgrid',false,'fracturelist',1:numel(fracplanes));
%-Frac-Matrix NNCs
G = fracturematrixShaleNNC3D(G,tol);
%-Frac-Frac NNCs
[G,fracplanes]=fracturefractureShaleNNCs3D(G,fracplanes,tol);
%-Well-Frac NNCs
[G,wells] = wellfractureNNCs3D(G,fracplanes,wells,tol);
G_edfm = G;

%% OMO: Projection-based NNCs
G = pMatFracNNCs3D(G,tol);
%%
% MRST includes both natural variables and overall composition. 
useNatural = true;

% Name of problem and pressure range
casename = 'c1c2c3';
pwf = 1000*psia;

%% Set up pEDFM operators
TPFAoperators = setupPEDFMOpsTPFA(G, G.rock, tol);

%% Define three-phase compressible flow model
nkr=2;
[fluid, info] = getShaleCompFluidCase(casename);
info.pressure = 5700*psia;  %this is redundant for the bakken case
info.temp = 352.594;

eosname = 'prcorr';  %'srk','rk','prcorr'
G1cell = cartGrid([1 1],[1 1]);
G1cell = computeGeometry(G1cell);
EOSModel = EquationOfStateModel(G1cell, fluid, eosname);
%Surface Conditions
p_sc = 101325; %atmospheric pressure
T_sc = 288.706;% 60 Farenheit
[~, ~, ~, ~, ~, rhoO_S, rhoG_S] = standaloneFlash(p_sc, T_sc, info.initial, EOSModel);
flowfluid = initSimpleADIFluid('n', [nkr, nkr, nkr], 'rho', [1000, rhoO_S, rhoG_S]);  

gravity reset off
model = NaturalVariablesCompositionalModel(G, [], flowfluid, fluid, 'water', false);
model.operators = TPFAoperators;

s0 = [1,0];
    
W = [];
for wi=1:numel(wells)
    W = addWellShale(W, G, G.rock, wells(wi).XFracCellIDs, ...
        'comp_i', [1,0],'Name', ['Prod',num2str(wi)], 'Val',...
        pwf, 'sign', -1, 'Type', 'bhp','Radius', wells(wi).radius,'Dir', 'x');
end
for wi=1:numel(wells)
    W(wi).components = info.initial;
end

plotWell(G,W);

M = csvread('CMG_timestep2.csv',1);    % Benchmark timesteps
dt_list=M(:,1)*day;
time_list=cumsum(convertTo(dt_list,day));

schedule = simpleSchedule(dt_list, 'W', W);

%% Impose initial pressure equilibrium
state = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, model.EOSModel);

%% Run simulations
[ws_pedfm , states_pedfm, report_pedfm] = simulateScheduleAD(state, model, schedule);
G_pedfm = G;






%% ------------------------------------------------------------------------
%                        EDFM Model Modifications
%  ------------------------------------------------------------------------
G=G_edfm;  % replace G with the previously-stored EDFM grid

%% Set up EDFM operators
TPFAoperators = setupShaleEDFMOpsTPFA(G, G.rock, tol);

%% Define three-phase compressible flow model
model = NaturalVariablesCompositionalModel(G, [], flowfluid, fluid, 'water', false);
model.operators = TPFAoperators;
    
W = [];
for wi=1:numel(wells)
    W = addWellShale(W, G, G.rock, wells(wi).XFracCellIDs, ...
        'comp_i', [1,0],'Name', ['Prod',num2str(wi)], 'Val',...
        pwf, 'sign', -1, 'Type', 'bhp','Radius', wells(wi).radius,'Dir', 'x');
end
for wi=1:numel(wells)
    W(wi).components = info.initial;
end

schedule = simpleSchedule(dt_list, 'W', W);

%% Impose initial pressure equilibrium
state = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, model.EOSModel);

%% Run simulations
[ws_edfm , states_edfm, report_edfm] = simulateScheduleAD(state, model, schedule);
G_edfm = G;







%% ------------------------------------------------------------------------
%                      Explicit Model Modifications
%  ------------------------------------------------------------------------
%% Create Box Reservoir Grid  
[G_matrix,frac_ids,well_ids]=explicitStencil3Fracs(physdim, 'aperture',w_f,...
      'numStages',1,'fracSpacing',physdim(1), 'fracHalfLength'  , x_f,      ...
       'nxRefine',12,'nxRefineSmall',10,'tipNY',8,'ny',8, 'nz',6, 'gridType','geomspacing');    

G_matrix = computeGeometry(G_matrix);

if (opt.shouldPlot)
    figure,
    plotGrid(G_matrix)%, view(5,45)
end

G_matrix.rock = makeRock(G_matrix, matrixPerm, matrixPoro);

%% Create Fracture System
G_matrix.rock.perm(frac_ids(1,:))=barrierPerm;
G_matrix.rock.perm(frac_ids(2,:))=k_f;
G_matrix.rock.perm(frac_ids(3,:))=barrierPerm;

%% Create Wells
wells = struct;
wells(1).XFracCellIDs=well_ids;
wells(1).radius=wellRadius;

%% visualize to check before pre-process
G=G_matrix;

%% Define shale gas flow model
gravity reset off
model = NaturalVariablesCompositionalModel(G, G.rock, flowfluid, fluid, 'water', false);

%% Assume constant BHP horizontal well
W = [];
for wi=1:numel(wells)
    W = addWell(W, G, G.rock, wells(wi).XFracCellIDs, ...
        'comp_i', [1,0],'Name', ['Prod',num2str(wi)], 'Val',...
        pwf, 'sign', -1, 'Type', 'bhp','Radius', wells(wi).radius,'Dir', 'x');
end

for i = 1:numel(W)
    W(i).components = info.initial;
end
schedule = simpleSchedule(dt_list, 'W', W);

%% Impose initial pressure equilibrium
state = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, model.EOSModel);

%% Run simulations
[ws, states, report] = simulateScheduleAD(state, model, schedule);

%% Plot the pEDFM, EDFM and Explicit Solution Profiles 

figure, 
plotToolbar(G_pedfm, states_pedfm)
view(40,30);
axis tight equal;

figure, 
plotToolbar(G_edfm, states_edfm)
view(40,30);
axis tight equal;

figure, 
plotToolbar(G, states)
view(40,30);
axis tight equal;

ws = {ws, ws_pedfm, ws_edfm};
states = {states, states_pedfm, states_edfm};

names = {'explicit', 'pEDFM', 'EDFM'};
plotWellSols(ws, cumsum(schedule.step.val), 'datasetnames', names,'field','qOs','linestyles',{'r','b'})

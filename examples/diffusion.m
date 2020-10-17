clear 
clc
close all
opt = struct('nkr',        1, ...
    'shouldPlot', 1 );

% fractureModel = 'EDFM'; 
fractureModel = 'pEDFM'; 

%% Load necessary modules, etc
mrstModule add hfm;             % hybrid fracture module
mrstModule add ad-core;
mrstModule add ad-props ad-core
mrstModule add mrst-gui;        % plotting routines
mrstModule add compositional; %mrstModule add ad-blackoil;
mrstModule add upr;

%% Basic Parameters
physdim = [1080 250 100]*meter;
griddim = [61 11 4];
cellsize = physdim./griddim;

%% Create Box Reservoir Grid
G_matrix = tensorGrid(0:cellsize(1):physdim(1),...
                      0:cellsize(2):physdim(2),...
                      0:cellsize(3):physdim(3));
G_matrix = computeGeometry(G_matrix);

if (opt.shouldPlot)
    figure,
    plotGrid(G_matrix), view(5,45)
end

%Rock Properties
matrix_perm = 1*nano*darcy;
matrix_poro= 0.063;
frac_aperture=0.0030;
frac_poro = 0.5;
frac_perm = 1*darcy;


G_matrix.rock=makeShaleRock(G_matrix, matrix_perm, matrix_poro);

fracCentroid1= [80.0, 125.0, 50.0]; % Centroid of the first fracture

%% Create Fracture System


[fracplanes,~] = createMultiStageHFs('numStages',10,'fracSpacing', 100,...
      'numFracsPerStage', 1,'fracHalfLength', 90,'fracHeight',60, ...
      'clusterSpacing', 10,'fracCentroid1',fracCentroid1,'isStencil',0);

for i=1:numel(fracplanes)
     fracplanes(i).aperture = frac_aperture;
     fracplanes(i).poro = frac_poro;
     fracplanes(i).perm = frac_perm;
end

if (opt.shouldPlot)
    figure,
    plotfracongrid(G_matrix,fracplanes,'label',false); % visualize to check before pre-process
end

%% Create Wells
wells = struct;
wells(1).points=[80.0 125.0, 50.0; 1000,125,50]*meter;
wells(1).radius=0.1*meter;

%% visualize to check before pre-process
G=G_matrix;
tol=1e-5;
%% EDFM PreProcessing

[G,fracplanes]=EDFMshalegrid(G,fracplanes,...
        'Tolerance',tol,'plotgrid',false,'fracturelist',1:numel(fracplanes));
  %-Frac-Matrix NNCs
G = fracturematrixShaleNNC3D(G,tol);
  %-Frac-Frac NNCs
[G,fracplanes]=fracturefractureShaleNNCs3D(G,fracplanes,tol);
  %-Well-Fracs NNCs
[G,wells] = wellfractureNNCs3D(G,fracplanes,wells,tol);

% OMO: Projection-based NNCs
if  strcmp(fractureModel,'pEDFM')
    G = pMatFracNNCs3D(G,tol);
end
if (opt.shouldPlot)
    figure,
    plotFracSystem(G,fracplanes,wells,'label',false)
end

% Name of problem and pressure range
casename = 'barnett3comps';
pwf = 500*psia;
pinj = 8000*psia;
G2=G;

%Use this to turn sorption on/off
% G.rock.shaleMechanisms.sorption = 1;

%diffusion
G.rock.shaleMechanisms.diffusion = 1;
G.rock.Di=[2.8,2.5,1.9]*10^-7;
G.rock.tau = 2;

G1=G;


%% Set up pEDFM operators
G=G1;
if  strcmp(fractureModel,'pEDFM')
    TPFAoperators = setupPEDFMOpsTPFA(G, G.rock, tol);
else
    TPFAoperators = setupEDFMOpsTPFA(G, G.rock, tol);
end


%% Define three-phase compositional flow model
% We define a three-phase compositional model with or without water.
% This is done by first instantiating the compositional model, and then
% manually passing in the internal transmissibilities and the topological
% neighborship from the embedded fracture grid.

[fluid, info] = getShaleCompFluidCase(casename);

eosname = 'prcorr';  %'srk','rk','prcorr'
G1cell = cartGrid([1 1],[1 1]);
G1cell = computeGeometry(G1cell);
EOSModel = EquationOfStateModel(G1cell, fluid, eosname);

%Surface Conditions
p_sc = 101325; %atmospheric pressure
T_sc = 288.706;% 60 Farenheit
[~, ~, ~, ~, ~, rhoO_S, rhoG_S] = standaloneFlash(p_sc, T_sc, info.initial, EOSModel);

flowfluid = initSimpleADIFluid('phases', 'WOG', 'n', [opt.nkr, opt.nkr, opt.nkr], 'rho', [1000, rhoO_S, rhoG_S]); 




gravity reset on

model = NatVarsShaleModel(G, [], flowfluid, fluid, 'water', true);
model.operators = TPFAoperators;

ncomp = fluid.getNumberOfComponents();
s0 = [0.23, 0.70, 0.07];

%% Set up initial state and schedule
% We set up a initial state with the reference pressure and a mixture of
% water, oil and gas initially. We also set up a simple-time step strategy 
% that ramps up gradually towards 30 day time-steps.



state1 = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, model.EOSModel);

W1 = [];
for wi=1:numel(wells)
    W1 = addWellShale(W1, G, G.Matrix.rock, wells(wi).XFracCellIDs, ...
        'comp_i', [0.23, 0.76, 0.01],'Name', ['Prod',num2str(wi)], 'Val',...
        pwf, 'sign', -1, 'Type', 'bhp','Radius', wells(wi).radius,'Dir','x');
end
for wi=1:numel(wells)
    W1(wi).components = info.initial;
end

totTime = 15*year;
nSteps = 15;
dt = rampupTimesteps(totTime, 30*day, nSteps);



schedule1 = simpleSchedule(dt, 'W', W1);

[ws1, states1, reports1] = simulateScheduleAD(state1, model, schedule1, 'Verbose', true);


%% Modeling system without sorption

clear model,
clear model.operators;
clear G;

G=G2;
if  strcmp(fractureModel,'pEDFM')
    TPFAoperators= setupPEDFMOpsTPFA(G, G.rock, tol);
else
    TPFAoperators= setupEDFMOpsTPFA(G, G.rock, tol);
end

% if useNatural
model = NaturalVariablesCompositionalModel(G, [], flowfluid, fluid, 'water', true);
% else
%     model = OverallCompositionCompositionalModel(G, [], flowfluid, fluid, 'water', true);
% end
model.operators = TPFAoperators;

state2 = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, model.EOSModel);

W2 = [];
for wi=1:numel(wells)
    W2 = addWellShale(W2, G, G.Matrix.rock, wells(wi).XFracCellIDs, ...
        'comp_i', [0.23, 0.76, 0.01],'Name', ['Prod',num2str(wi)], 'Val',...
        pwf, 'sign', -1, 'Type', 'bhp','Radius', wells(wi).radius,'Dir','x');
end
for wi=1:numel(wells)
    W2(wi).components = info.initial;
end

schedule2= simpleSchedule(dt, 'W', W2);

[ws2, states2, reports2] = simulateScheduleAD(state2, model, schedule2, 'Verbose', true);

%% plotting

figure, 
plotToolbar(G1, states1)
view(40,30);
axis tight equal;

figure, 
plotToolbar(G2, states2)
view(40,30);
axis tight equal;

ws = {ws1, ws2};
states = {states1, states2};

names = {'With Difussion', 'Without Difussion'};
plotWellSols(ws, cumsum(schedule2.step.val), 'datasetnames', names,'field','qGs','linestyles',{'r','b'})

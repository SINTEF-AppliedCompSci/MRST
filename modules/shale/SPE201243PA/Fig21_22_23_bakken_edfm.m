%% This will generate Figures 21-23 in SPE-201243-PA
%
clear 
clc
close all
Globals; 
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
physdim = [1200 300 90]*meter;
griddim = [40 12 5];
cellsize = physdim./griddim;
theta=90;
theta_1=theta*(pi/180);

%% Create Fracture System
[fracplanes,frac_centroid_s]=MulistagePlanarNonPlanarHF('numStages',8,'fracSpacing', 100,...
      'numFracsPerStage', 3,'fracHalfLength', 70,'fracHeight', 60,'theta',theta_1,...
      'clusterSpacing', 10,'heelCoord', [125.0 150.0 45.0]);
  
for i=1:numel(fracplanes)
     fracplanes(i).aperture = 3.0*milli*meter; %1*micro*meter;
     fracplanes(i).poro = 0.5;
     fracplanes(i).perm = 1.0*darcy;%100*micro*darcy;
end

G_matrix = StructuredGridFractureonFace(physdim,frac_centroid_s,'numStages',8,'numFracsPerStage', 3,...
    'ny',griddim(2),'nz',griddim(3),'nxRefine_small',4,'fracSpacing', 100,'nxRefine',5,'aperture', 3.0*milli*meter);

G_matrix = computeGeometry(G_matrix);


G_matrix.rock=makeRock(G_matrix,100*nano*darcy,0.063);
if (opt.shouldPlot)
    plotfracongrid(G_matrix,fracplanes,'label',false); % visualize to check before pre-process
end
%% Stochastic Natural Fracture Generation using ADFNE
% Globals; 
rng(12345); 
tol=1e-5;
% tol4domain = max(cellsize)+tol;%This is to prevent having fractures on the boundary
tol4domain = tol;

set1 = Field(DFN('dim',3,'n',500,'dir',45,'ddir',-1e9,'minl',10,...
            'mu',10,'maxl',20,'bbx',[tol4domain,tol4domain,tol4domain,physdim(1)-tol4domain,physdim(2)-tol4domain,physdim(3)-tol4domain],'dip',45,'ddip',-1e9,...
            'shape','l','q',4),'Poly');  
set2 = Field(DFN('dim',3,'n',500,'dir',-45,'ddir',-1e9,'minl',10,...
            'mu',10,'maxl',20,'bbx',[tol4domain,tol4domain,tol4domain,physdim(1)-tol4domain,physdim(2)-tol4domain,physdim(3)-tol4domain],'dip',45,'ddip',-1e9,...
            'shape','l','q',4),'Poly');


[set1_,nonPlanarSets1,fracArea1] = processStochFracs(set1);
[set2_,nonPlanarSets2,fracArea2] = processStochFracs(set2); 
fracArea = [fracArea1;fracArea2];

fprintf('%d of %d set1 fracs were OK while %d of %d set2 fracs were OK \n',...
    numel(set1_),numel(set1),numel(set2_),numel(set2));

fprintf('Number of nonplanar sets in sets 1 and 2 are : %d and %d respectively\n',...
    numel(nonPlanarSets1),numel(nonPlanarSets2));

figure;
if (opt.shouldPlot)
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
    fracplanes(idxGlobal).aperture = 100*nano*meter; %1*micro*meter;
    fracplanes(idxGlobal).poro=0.5;
    if(mod(i,2)==false)%even NF idx
       fracplanes(idxGlobal).perm=100*micro*darcy;%100*micro*darcy;%0.01*darcy;
    else%odd NF idx
       fracplanes(idxGlobal).perm=1*nano*darcy;%100*micro*darcy;%0.01*darcy;
    end
end

% if (opt.shouldPlot)
%     figure,
%     plotfracongrid(G_matrix,fracplanes,'label',false); % visualize to check before pre-process
% end

%% Create Wells
wells = struct;
wells(1).points=[90.0 150.0 45.0; 1000,150,45]*meter;
wells(1).radius=0.1*meter;

%% visualize to check before pre-process
G=G_matrix;

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
% G = pMatFracNNCs3D(G,tol);


if (opt.shouldPlot)
    figure,
    plotfracSystem_NF(G,fracplanes,numHFplanes,wells,'label',false)
end

%%
useNatural = true;

% Name of problem and pressure range
casename = 'bakken';
pwf = 2000*psia;
pinj = 8000*psia;


%% Set up EDFM operators
TPFAoperators = setupShaleEDFMOpsTPFA(G, G.rock, tol);

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

if useNatural
    model = NaturalVariablesCompositionalModel(G, [], flowfluid, fluid, 'water', true);
else
    model = OverallCompositionCompositionalModel(G, [], flowfluid, fluid, 'water', true);
end
model.operators = TPFAoperators;

%% Set up initial state and schedule
% We set up a initial state with the reference pressure and a mixture of
% water, oil and gas initially. We also set up a simple-time step strategy 
% that ramps up gradually towards 30 day time-steps.

% totTime = 15*year;
% nSteps = 25;

ncomp = fluid.getNumberOfComponents();
s0 = [0.23, 0.70, 0.07];   

state = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, model.EOSModel);

W = [];
for wi=1:numel(wells)
    W = addWellShale(W, G, G.Matrix.rock, wells(wi).XFracCellIDs, ...
        'comp_i', [0.23, 0.76, 0.01],'Name', ['Prod',num2str(wi)], 'Val',...
        pwf, 'sign', -1, 'Type', 'bhp','Radius', wells(wi).radius,'Dir','x');
end
for wi=1:numel(wells)
    W(wi).components = info.initial;
end
% plotWell(G,W);
% M = csvread('CMG_timestep2.csv',1);%Benchmark_CMG\CMG_timestep.csv
% dt_list=M(:,1)*day;
% schedule = simpleSchedule(dt_list, 'W', W);
% time1=dt_list/86400;
% time=cumsum(time1);%for plotting in excel 
% schedule = simpleSchedule(dt, 'W', W);
totTime = 30*year;
nSteps = 45;
dt = rampupTimesteps(totTime, 1*year, nSteps);
schedule = simpleSchedule(dt, 'W', W);

time=cumsum(dt)/86400;%for plotting in excel 


%% Simulate problem
[ws, states, reports] = simulateScheduleAD(state, model, schedule, 'Verbose', true);

%% plotting

figure, 
plotToolbar(G, states)
view(40,30);
axis tight equal;

figure,
plotWellSols(ws,cumsum(schedule.step.val))

%Export it to ParaviewResults Folder
% Write2Paraview([mfilename '.vtk'],G,fracplanes,wells,states);
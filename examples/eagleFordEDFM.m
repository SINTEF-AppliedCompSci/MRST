clc;
clear;
close all;
try
    %% This code requires ADFNE and the Statistics & Machine Learning Toolbox
    % The ADFNE folder should be put within the shale module folder
    Globals;
    warning('on');
catch
    %% Print warning when ADFNE is not available
    warning('adfne:missing', 'ADFNE or Statistics & Machine Learning Toolbox is not available!');
    warning('off', 'adfne:missing');
end

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
tol=1e-5;
physdim = [150 300 38.1];
fracHalfLength=180;
fracHeight=75*ft;
fractureSpacing=150;
fracPerStage=3;
NumStages=1;
clusterSpacing=10;
fracAperture=3.0*milli*meter;
fracPoro=0.5;
fracPerm=1.0*darcy;
MatrixPoro=0.12;
MatrixPerm=900*nano*darcy;
wellLength=150;
wellRadius=0.328084*ft;
heelCooord=[0, 0, physdim(3)-tol];
ToeCoord=[heelCooord(1)+wellLength,0, physdim(3)-tol];
BHP=1000*psia;
fracCentroid1= [65, fracHalfLength/2, physdim(3)-fracHeight/2]; % Centroid of the first fracture
%% Create Fracture System
[fracplanes,frac_centroid_s] = createMultiStageHFs('numStages',NumStages,'fracSpacing', fractureSpacing,...
      'numFracsPerStage', fracPerStage,'fracHalfLength', fracHalfLength,'fracHeight',fracHeight, ...
      'clusterSpacing', clusterSpacing,'fracCentroid1',fracCentroid1,'isStencil',1);
  
for i=1:numel(fracplanes)
     fracplanes(i).aperture = fracAperture; 
     fracplanes(i).poro = fracPoro;
     fracplanes(i).perm = fracPerm;
end
G_matrix = meshHFsystem(physdim,frac_centroid_s,'numStages',NumStages,'numFracsPerStage',fracPerStage,...
    'ny',20,'nz',5,'nxRefine_small',4,'nxRefine',9,'fracSpacing', fractureSpacing,'aperture',fracAperture);


G_matrix = computeGeometry(G_matrix);

G_matrix.rock=makeRock(G_matrix,MatrixPerm,MatrixPoro);
if (opt.shouldPlot)
    plotfracongrid(G_matrix,fracplanes,'label',false); % visualize to check before pre-process
end
view(40,30);


%% Stochastic Natural Fracture Generation using ADFNE
rng(12345); 


try
    %% This code requires ADFNE and the Statistics & Machine Learning Toolbox
    set1 = Field(DFN('dim',3,'n',75,'dir',45,'ddir',-100,'minl',5,'mu',20, ...
        'maxl',20,'bbx',[tol,tol,tol,physdim(1)-tol, ...
        physdim(2)-tol,physdim(3)-tol],'dip',45,'ddip',-100, ...
        'shape','l','q',4),'Poly');
    set2 = Field(DFN('dim',3,'n',75,'dir',-45,'ddir',-100,'minl',5,'mu',20, ...
        'maxl',20,'bbx',[tol,tol,tol,physdim(1)-tol, ...
        physdim(2)-tol,physdim(3)-tol],'dip',45,'ddip',-100, ...
        'shape','l','q',4),'Poly');
catch
    %% Load generated fractures when ADFNE is not available
    load('EagleFordFracs.mat');
end
    
[set1_,nonPlanarSets1,fracArea1] = processStochFracs(set1);
[set2_,nonPlanarSets2,fracArea2] = processStochFracs(set2); 
fracArea = [fracArea1;fracArea2];

fprintf('%d of %d set1 fracs were OK while %d of %d set2 fracs were OK \n',...
    numel(set1_),numel(set1),numel(set2_),numel(set2));

fprintf('Number of nonplanar sets in sets 1 and 2 are : %d and %d respectively\n',...
    numel(nonPlanarSets1),numel(nonPlanarSets2));

figure; 
if (opt.shouldPlot) %ADFNE Plot funcs     
    try      
    % Plot stereonet if ADFNE is available
        figure;          
        Draw('ply',set1_);    
        Draw('ply',set2_);view(45,30)             
        o = Orientation([set1_;set2_]);          
        figure,          
        Stereonet([o.Dip],[o.Dir],'density',true,'marker','*','ndip',6,...              
            'ndir',24,'cmap',@jet,'color','y');   
    catch 
    % Don’t plot if ADFNE is not available             
    end    
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
G_matrix.numHFplanes=numHFplanes;
G_matrix.numNFplanes=numNFplanes;

%% Create Wells
wells = struct;
wells(1).points=[heelCooord; ToeCoord];
wells(1).radius=wellRadius;

%% visualize to check before pre-process
G=G_matrix;

%% EDFM PreProcessing
% The preprocessing takes a while. The user can  remove line 147 to 157 and
% instead load the preprocessed grid at line 160

[G,fracplanes]=EDFMshalegrid(G,fracplanes,...
        'Tolerance',tol,'plotgrid',false,'fracturelist',1:numel(fracplanes));
  %-Frac-Matrix NNCs
G = fracturematrixShaleNNC3D(G,tol);
  %-Frac-Frac NNCs
[G,fracplanes]=fracturefractureShaleNNCs3D(G,fracplanes,tol);
  %-Well-Fracs NNCs
[G,wells] = wellfractureNNCs3D(G,fracplanes,wells,tol);

% Set up EDFM operators
TPFAoperators = setupShaleEDFMOpsTPFA(G, G.rock, tol);
 
% Load preprocessed grid
% load('EagleFordEDFMGrid.mat');

if (opt.shouldPlot)
    figure,
    plotfracSystem_NF(G,fracplanes,numHFplanes,wells,'label',false)
end
view(40,30);

% Name of problem and pressure range
casename = 'eagleford';
pwf = BHP;


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

%% Set up initial state and schedule

ncomp = fluid.getNumberOfComponents();
s0 = [0.17, 0.74, 0.09];   

state = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, model.EOSModel);

W = [];
for wi=1:numel(wells)
    W = addWellShale(W, G, G.Matrix.rock, wells(wi).XFracCellIDs, ...
        'comp_i', [0.17, 0.74, 0.09],'Name', ['Prod',num2str(wi)], 'Val',...
        pwf, 'sign', -1, 'Type', 'bhp','Radius', wells(wi).radius,'Dir','x');
end
for wi=1:numel(wells)
    W(wi).components = info.initial;
end

totTime = 15*year;
nSteps = 15;
dt = rampupTimesteps(totTime, 30*day, nSteps);

schedule = simpleSchedule(dt, 'W', W);
%% Simulate problem
[ws, states, reports] = simulateScheduleAD(state, model, schedule, 'Verbose', true);

%% plotting

figure, 
plotToolbar(G, states)
view(40,30);
axis tight equal;

figure,
plotWellSols(ws,cumsum(schedule.step.val))

%% Loading Modules
clear; 
clc;
close all;
Globals;
case2run ='ProdBot';%'ProdBot_InjTop'; %'ProdBot';
opt = struct('nkr',        1, ...
             'shouldPlot', 0 ); %change to 0 if running on HPC
 %% Load necessary modules, etc 
mrstModule add hfm;             % hybrid fracture module
mrstModule add ad-core;         
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
dz=6000*ft*ones([numel(0:5:physdim(1)), numel(0:5:physdim(2))]);
G_matrix = tensorGrid(0:5:physdim(1), 0:5:physdim(2),z, 'depthz', dz);
G_matrix = computeGeometry(G_matrix);    
G_matrix.rock=makeRock(G_matrix,10*micro*darcy,0.07); %20*micro*darcy
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
[fracIndx_X,fracIndx_Y,fracIndx_Z] = meshgrid(1:G_matrix.cartDims(1), 1:G_matrix.cartDims(2), frac_z);
fraccells = sub2ind(G_matrix.cartDims, reshape(fracIndx_X,numel(fracIndx_X),1),reshape(fracIndx_Y,...
            numel(fracIndx_X),1), reshape(fracIndx_Z,numel(fracIndx_X),1));
G_matrix.rock.poro(fraccells) = 0.33;% 0.33*3/100;  %should be related to size of w_f wrt frac cell size
G_matrix.rock.perm(fraccells) = 10*darcy; %5*darcy

if (opt.shouldPlot)
    figure,
    plotCellData(G_matrix,convertTo(G_matrix.rock.perm,milli*darcy));
    colorbar('horiz'); axis equal tight; view(3);
end

G=G_matrix;
%% Define three-phase compressible flow model
useNatural = true;
casename = 'oil_1'; %'bakken_light';
pwf = 2500*psia;
pinj = 3000*psia;
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
%     nkr = 1;
flowfluid = initSimpleADIFluid('phases', 'OG', 'n', [opt.nkr, opt.nkr], 'rho', [ rhoO_S, rhoG_S]);    % flowfluid.KGangiFn = @(p) power((1-power(((Pc - alpha.*p)./Pmax),m)),3);
gravity reset on

arg = {G, G.rock, flowfluid, fluid, 'water', false};

diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);
mex_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'useMex', true);
sparse_backend = SparseAutoDiffBackend();

if useNatural
    constructor = @NaturalVariablesCompositionalModel;
else
    constructor = @GenericOverallCompositionModel;
end

modelSparseAD = constructor(arg{:}, 'AutoDiffBackend', sparse_backend);
modelDiagonalAD = constructor(arg{:}, 'AutoDiffBackend', diagonal_backend);
modelMexDiagonalAD = constructor(arg{:}, 'AutoDiffBackend', mex_backend);

%% Set up initial state
totTime = 4*year;%60*year
nSteps =15;
ncomp = fluid.getNumberOfComponents();
s0 = [ 1, 0];   %s0 = [0.23, 0.77, 0.07];
%                                 (G, p, T, s0, z0, eos)
state = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, modelSparseAD.EOSModel);
%% Pick linear solver
% The AMGCL library is one possible solver option for MRST. It is fairly
% easy to write interfaces to other solvers using MEX files and/or the
% LinearSolverAD class. We call the routine for automatically selecting a
% reasonable linear solver for the specific model.
linsolve = selectLinearSolverAD(modelSparseAD,'useAMGCL',true,'useCPR',true);
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
    case 'ProdBot' 
        % Producer
        wellCell = sub2ind(G_matrix.cartDims,20,10,frac_z(2));
        W = addWell(W, G, G.rock, wellCell, ...
            'comp_i', [1, 0],'Name', 'Prod_Bot', 'Val', pwf, 'sign', -1,...
            'Type', 'bhp','Radius', wellRadius,'Dir','z','refDepth', 6000*ft+70);
        W(1).components = info.initial;
    otherwise
        warning('Case Does Not Exist. Running case with Prod Only at Bottom')
        % Producer
        W = verticalWell(W, G.Matrix, G.Matrix.rock, 1, 1, 80, ...
            'comp_i', [0.23, 0.76, 0.01],'Name', 'Prod_Bot', 'Val', pwf, 'sign', -1, 'Type', 'bhp','Radius', wellRadius);
        W(1).components = info.initial;
end
plotWell(G,W); 
dt = rampupTimesteps(totTime, 30*day, nSteps); %20*day
schedule = simpleSchedule(dt, 'W', W);  
%% Simulate problem
[ws, states, reports] = simulateScheduleAD(state, modelMexDiagonalAD, schedule, 'nonlinearsolver', nls, 'Verbose', true);
%% plotting
figure, 
plotToolbar(G, states)
view(40,30);
axis tight equal;
plotWellSols(ws,cumsum(schedule.step.val))
%% Calculate RF
tinSecs = cumsum(schedule.step.val);
tinDays = tinSecs./86400;
tinYears = tinDays./365.25;

[wsEDFMflashed,qProd] = flashWS(ws, 'recoveryMethod', 'SDEOR','casename', 'oil_1');
Nprod = cumtrapz(tinSecs,qProd);
plotWellSols(wsEDFMflashed,cumsum(schedule.step.val))

pWellTime = zeros(length(ws),1);
for i=1:length(ws)
    pAllCells = states{i}.pressure/psia;
    pWellTime(i) = pAllCells(wellCell);
end

InitialMass = sum(states{1}.rho(:,1).*G.cells.volumes.*G.rock.poro);
CumMassProd1 = zeros(64,1);
for i=1:64
    CumMassProd1(i)=InitialMass-sum(states{i}.rho(:,1).*G.cells.volumes.*G.rock.poro);%total mass produced in Kg
end
CumMassProd2 = cumtrapz(tinSecs,qProd); %total mass produced in Kg
figure,
plot(tinYears,CumMassProd1,tinYears,-CumMassProd2)
legend('CumMassProd','CumMassProd2')
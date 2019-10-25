%% Upscaling workflow
mrstModule add ad-core ad-blackoil ad-props mrst-gui msrsb
close all;

%% Build simulation model
% Define areal mesh
[xmax,ymax, n]  = deal(1000*meter, 1000*meter, 30);
[x, y] = meshgrid(linspace(0,xmax,n+1), linspace(0,ymax,n+1));
[x, y] = deal(x',y');

% Basic dome structure
dome = 1-exp(sqrt((x - xmax/2).^2 + (y - ymax/2).^2)*1e-3);

% Periodic perturbation
[xn,yn] = deal(pi*x/xmax,pi*y/ymax);
perturb = sin(5*xn) + .5*sin(4*xn+6*yn) + cos(.25*xn*yn./pi^2) + cos(3*yn);
perturb = perturb/3.5;

% Random small-scale perturbation
rng(0);
[h, hr] = deal(8,1);
zt = 50 + h*perturb + rand(size(x))*hr - 20*dome;
zb = zt + 30;
zmb = min(zb + 4 + 0.01*x - 0.020*y + hr*rand(size(x)), zb);
zmt = max(zb -15 + 0.01*x - 0.025*y + hr*rand(size(x)), zt);
horizons = {struct('x', x, 'y', y, 'z', zt), ...
    struct('x', x, 'y', y, 'z', zmt), ...
    struct('x', x, 'y', y, 'z', zmb), ...
    struct('x', x, 'y', y, 'z', zb)};


% Interpolate to build unfaulted corner-point grid
dims = [40, 40]; layers = [3 6 3];
grdecl = convertHorizonsToGrid(horizons, 'dims', dims, 'layers', layers);

% Insert faults
[X,Y,Z]  = buildCornerPtNodes(grdecl);

i=47:80; Z(i,:,:) = Z(i,:,:) + .022*min(0,Y(i,:,:)-550);
j= 1:30; Z(:,j,:) = Z(:,j,:) + .021*min(0,X(:,j,:)-400);
j=57:80; Z(:,j,:) = Z(:,j,:) + .023*min(0,X(:,j,:)-750);
grdecl.ZCORN = Z(:);

G = processGRDECL(grdecl);
G = computeGeometry(G);

% Petrophysics
% Set up permeability based on K-indices and introduce anisotropy by
% setting K_z = .1*K_x
rng(357371);
[K,L] = logNormLayers(G.cartDims, [100 400 10 50]*milli*darcy);
K = K(G.cells.indexMap);
perm = [K, K, 0.1*K];
rock = makeRock(G, perm, 0.3);

% Define wells
% Producers
simTime = 10*year;
pv      = poreVolume(G, rock);
injRate = 1*sum(pv)/simTime;
offset  = 10;
W = verticalWell([], G, rock, offset, offset, [],...
                'Name', 'P1', 'comp_i', [1 0], ...
                'Val', 250*barsa, 'Type', 'bhp', 'refDepth', 50);
W = verticalWell(W, G, rock,  offset, floor(G.cartDims(1)/2)+3, [],...
                'Name', 'P2', 'comp_i', [1 0], ...
                'Val', 250*barsa, 'Type', 'bhp', 'refDepth', 50);
W = verticalWell(W, G, rock, offset, G.cartDims(2) - offset/2, [], ...
                'Name', 'P3', 'comp_i', [1 0], ...
                'Val', 250*barsa, 'Type', 'bhp', 'refDepth', 50);

% Injectors
W = verticalWell(W, G, rock, G.cartDims(1)-5, offset, [],...
                'Name', 'I1', 'comp_i', [1 0], ...
                'Val', injRate, 'Type', 'rate', 'refDepth', 50);


% Define fluid behavior and instantiate model
% Three-phase template model
fluid = initSimpleADIFluid('mu',    [1, 5, 0]*centi*poise, ...
                           'rho',   [1000, 700, 0]*kilogram/meter^3, ...
                           'n',     [2, 2, 0]);

% Constant oil compressibility
c        = 0.001/barsa;
p_ref    = 300*barsa;
fluid.bO = @(p, varargin) exp((p - p_ref)*c);

% Construct reservoir model
gravity reset on
model = TwoPhaseOilWaterModel(G, rock, fluid);


% Define initial state
depthOW = 85*meter;
depthD  = 10*meter;
region = getInitializationRegionsBlackOil(model, depthOW, ...
            'datum_depth', depthD, 'datum_pressure', p_ref);
state0 = initStateBlackOilAD(model, region);


%% Define simulation schedule and set solver parameters

% Compute the timestep
nstep   = 25;
refine  = 5;
startSteps = repmat((simTime/(nstep + 1))/refine, refine, 1);
restSteps=  repmat(simTime/(nstep + 1), nstep, 1);
timesteps = [startSteps; restSteps];

% Set up the schedule containing both the wells and the timestep
schedule = simpleSchedule(timesteps, 'W', W);

% Tighten tolerences
model.drsMaxRel = inf;
model.dpMaxRel  = .1;
model.dsMaxAbs  = .1;

% Set up CPR preconditioner
try
    mrstModule add agmg
    pressureSolver = AGMGSolverAD('tolerance', 1e-4);
catch
    pressureSolver = BackslashSolverAD();
end
linsolve = CPRSolverAD('ellipticSolver', pressureSolver, 'relativeTolerance', 1e-3);


%% Simulate base case
fn = getPlotAfterStep(state0, model, schedule,'view',[50 50], ...
                     'field','s:1','wells',W);
[wellSols, states, report] = ...
   simulateScheduleAD(state0, model, schedule, ...
                       'LinearSolver', linsolve, 'afterStepFn',fn);

%% Simulate with multiscale solver
cdims  = [10, 10, 3];
p0     = partitionUI(G, cdims);
Gf     = makeInternalBoundary(G, find(G.faces.tag > 0));
p      = processPartition(Gf, p0);
CG     = coarsenGeometry(generateCoarseGrid(G, p));
CG     = storeInteractionRegion(CG, 'edgeBoundaryCenters', true, 'adjustCenters', true);
msSolver = MultiscaleVolumeSolverAD(CG, 'tolerance', 1e-6, ...
    'maxIterations', 100, 'useGMRES', true, ...
    'getSmoother', getSmootherFunction('type', 'ilu0', 'iterations', 1));

linsolve = CPRSolverAD('ellipticSolver', msSolver, 'relativeTolerance', 1e-3);

fn = getPlotAfterStep(state0, model, schedule,'view',[50 50], ...
                     'field','s:1','wells',W);
[wellSolsMS, statesMS, reportMS] = ...
   simulateScheduleAD(state0, model, schedule, ...
                       'LinearSolver', linsolve, 'afterStepFn',fn);

%% Plot comparison
plotWellSols({wellSols, wellSolsMS}, cumsum(schedule.step.val),'datasetnames',{'agmg','msrsb'});
%% Load modules
mrstModule add ad-core ad-blackoil ad-props deckformat mrst-gui linearsolvers
%% 
gravity reset on
mrstVerbose on
useMex = true;
opm = mrstPath('opm-tests');
assert(~isempty(opm), 'You must register https://github.com/opm/opm-tests as a module!');
[deck, output] = getDeckOPMData('norne', 'NORNE_ATW2013');
%% Build model from EGRID/INIT files
egrid = readEclipseOutputFileUnFmt([output.opm.location, '.EGRID']);
init = readEclipseOutputFileUnFmt([output.opm.location, '.INIT']);
[Ge, rock_ecl, Ne, Te] = initGridFromEclipseOutput(init, egrid, 'outputSimGrid', true);
G_viz = computeGeometry(Ge{1});
G_sim = Ge{2};


rock = initEclipseRock(deck);
rock = compressRock(rock, G_sim.cells.indexMap);

fluid = initDeckADIFluid(deck, 'useMex', useMex);
% Setup model, but skip setting up the operators since we do not have a
% proper grid
model = GenericBlackOilModel(G_sim, [], fluid, 'disgas', true, 'vapoil', true, 'inputdata', deck);
% Finally set up the connectivity graph from output
model.rock = rock;
model.operators = setupOperatorsTPFA(G_sim, rock, 'deck', deck, 'neighbors', Ne, 'trans', Te);
%% Set up everything
[state0, model, schedule, nls] = initEclipseProblemAD(deck, 'model', model, ...
                                                        'TimestepStrategy', ...
                                                        'ds', 'useCPR', true, ...
                                                        'useMex', useMex);
nls.LinearSolver.tolerance = 1e-3;
%% Pack problem
% Tweaks to have somewhat similar criterion as OPM. For my comparison I ran
% OPM with:
% ToleranceCnv="0.01" # default: "0.01"
% ToleranceCnvRelaxed="0.01" # default: "1e+09"
% ToleranceMb="1e-7" # default: "1e-06"

nls.maxIterations = 12;
% nls.timeStepSelector = SimpleTimeStepSelector();
nls.timeStepSelector = FactorTimeStepSelector();

model.toleranceCNV = 1e-2;
model.toleranceMB = 1e-7;

% Set well tolerances
model.FacilityModel = ExtendedFacilityModel(model);
model.FacilityModel.toleranceWellBHP = 1e-3;
model.FacilityModel.toleranceWellRate = 5e-3;

% Reset just in case
model.FluxDiscretization = [];
model.FlowPropertyFunctions = [];
model.PVTPropertyFunctions = [];

model.AutoDiffBackend.useMex = true;
model.AutoDiffBackend.rowMajor = true;

model = model.validateModel();
% Use the alternative more rigorous crossflow definition for component
% fluxes
xflow = WellComponentTotalVolumeBalanceCrossflow(model);
xflow.onlyLocalDerivatives = true;

model.FacilityModel.FacilityFluxDiscretization.ComponentTotalFlux = xflow;


useFlag = false;
model.PVTPropertyFunctions.Viscosity.useSaturatedFlag = useFlag;
model.PVTPropertyFunctions.ShrinkageFactors.useSaturatedFlag = useFlag;

RvMax = model.PVTPropertyFunctions.getStateFunction('RvMax');
RvMax.rvReduction = 0.5; % Not tested, but in the deck
model.PVTPropertyFunctions = model.PVTPropertyFunctions.setStateFunction('RvMax', RvMax);


model.dsMaxAbs = 0.2;
model.dpMaxRel = 0.1;

nls.useRelaxation = true;
nls.minRelaxation = 0.5;
nls.maxIterations = 18;
nls.maxTimestepCuts = 10;
nls.oscillationThreshold = 0.5;
nls.acceptanceFactor = 1;
nls.verbose = false;
nls.LinearSolver.verbose = false;

problem = packSimulationProblem(state0, model, schedule, 'norne', ...
    'Name', 'GenericBlackOil_FI_EGRID', 'nonlinearsolver', nls);

%% MRST grid     
simulatePackedProblem(problem, 'continueOnError', false);

%% Get results
[ws, states, reports] = getPackedSimulatorOutput(problem);
nstep = numel(ws);
T = cumsum(schedule.step.val(1:nstep));
%% Plot well sols

% Get output stored in test repository
ws_opm = output.opm.wellSols(2:end);
T_opm = cumsum(schedule.step.val);
fname = 'OPM-Flow-Legacy';
wellSols = {ws, ws_opm};
time = {T; T_opm};
names = {'MRST', fname};

try
    % Flow output - may not be present
    pp = fullfile(mrstPath('opm-tests'), 'norne', 'NORNE_ATW2013');
    states_newopm = convertRestartToStates(pp, model.G);
    [ws_opm, T_opm_new] = convertSummaryToWellSols(pp);
    fname = 'OPM-Flow';
    
    wellSols{end+1} = ws_opm;
    time{end+1} = T_opm_new;
    names{end+1} = fname;
catch
    
end

ecl = fullfile(mrstPath('opm-tests'), 'norne', 'ECL.2014.2', 'NORNE_ATW2013');
[ws_ecl, T_ecl] = convertSummaryToWellSols(ecl, 'metric');
wellSols{end+1} = ws_ecl;
time{end+1} = T_ecl;
names{end+1} = 'Eclipse';

try
    % Flow output without hysteresis to match MRST - may not be present
    fn = fullfile(mrstPath('scratchpad'), 'new_fluid', 'norne');
    rstrt = load(fullfile(fn, 'restart_nohyst.mat'));
    rstrt = rstrt.rs;
    % 
    smry = load(fullfile(fn, 'smry_nohyst'));
    smry = smry.smry;
    [ws_nohyst, T_nohyst] = convertSummaryToWellSols(smry, 'metric');
    wellSols{end+1} = ws_nohyst;
    time{end+1} = T_nohyst;
    names{end+1} = 'Flow-NoHyst';
catch
    
end
plotWellSols(wellSols, time, 'datasetnames', names, 'linestyles', {'-o', '--'})
%% Plot states
figure;
plotToolbar(G_viz, output.opm.states);
title('OPM');
figure;
plotToolbar(G_viz, states);
title('MRST');


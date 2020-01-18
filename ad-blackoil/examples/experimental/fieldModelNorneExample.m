%% Load modules
mrstModule add ad-core ad-blackoil ad-props deckformat mrst-gui linearsolvers
%% 
mrstVerbose on
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

fluid = initDeckADIFluid(deck);
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
                                                        'useMex', true);
%% Pack problem
% Tweaks to have somewhat similar criterion as OPM. For my comparison I ran
% OPM with:
% ToleranceCnv="0.01" # default: "0.01"
% ToleranceCnvRelaxed="0.01" # default: "1e+09"
% ToleranceMb="1e-7" # default: "1e-06"

nls.maxIterations = 12;
model.toleranceCNV = 1e-2;
model.toleranceMB = 1e-7;

% Set well tolerances
model.FacilityModel = ExtendedFacilityModel(model);
model.FacilityModel.toleranceWellBHP = 1e-3;
model.FacilityModel.toleranceWellRate = 1e-3;
% Use the alternative more rigorous crossflow definition for component
% fluxes
model.FacilityModel.FacilityFluxDiscretization.ComponentTotalFlux = WellComponentTotalFluxDensityMix(model);

% Reset just in case
model.FluxDiscretization = [];
model.FlowPropertyFunctions = [];
model = model.validateModel();
useFlag = false;
model.FlowPropertyFunctions.Viscosity.useSaturatedFlag = useFlag;
model.FlowPropertyFunctions.ShrinkageFactors.useSaturatedFlag = useFlag;

RvMax = model.FlowPropertyFunctions.getStateFunction('RvMax');
RvMax.rvReduction = 0.5; % Not tested, but in the deck
model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('RvMax', RvMax);


model.dsMaxAbs = 0.2;
model.dpMaxRel = 0.1;

nls.useRelaxation = true;
nls.minRelaxation = 0.5;

problem = packSimulationProblem(state0, model, schedule, 'norne', ...
    'Name', 'GenericBlackOil_FI_EGRID', 'nonlinearsolver', nls);

%% MRST grid     
simulatePackedProblem(problem, 'continueOnError', false);

%% Get results
[ws, states, reports] = getPackedSimulatorOutput(problem);
nstep = numel(ws);
T = cumsum(schedule.step.val(1:nstep));
%% Plot well sols
if 1
    % Get output stored in test repository
    ws_opm = output.opm.wellSols(2:end);
    T_opm = cumsum(schedule.step.val);
    fname = 'OPM-Flow-Legacy';
else
    pp = fullfile(mrstPath('opm-tests'), 'norne', 'NORNE_ATW2013');
    [ws_opm, T_opm] = convertSummaryToWellSols(pp);
    fname = 'OPM-Flow';
end
wellSols = {ws, ws_opm};

time = {T; T_opm};
names = {'MRST', fname};
plotWellSols(wellSols, time, 'datasetnames', names)
%% Plot states
figure;
plotToolbar(G_viz, output.opm.states);
title('OPM');
figure;
plotToolbar(G_viz, states);
title('MRST');


%% Example demonstrating how to run the Norne field model in MRST
% Note: The example uses a custom well function, and MRST does not support
% hysteresis. The model is read in as an EGRID and is initialized in MRST
%
% Norne is a real field model, released by Equinor and its partners. It is
% a fairly complex blackoil model with rel. perm. scaling, threshold
% pressures, many wells, historical resv controls and many other features
% not normally seen in academic simulations.
mrstModule add ad-core ad-blackoil ad-props deckformat mrst-gui linearsolvers
%% Set up case
gravity reset on
mrstVerbose on
useMex = true;
opm = mrstPath('opm-tests');
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
model.operators = setupOperatorsTPFA(G_sim, rock, 'neighbors', Ne, 'trans', Te);
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
nls.timeStepSelector = FactorTimeStepSelector();
% Same tolerances as OPM
model.toleranceCNV = 1e-2;
model.toleranceMB = 1e-7;

% Set well tolerances
model.FacilityModel = GenericFacilityModel(model);
model.FacilityModel.toleranceWellBHP = 1e-3;
model.FacilityModel.toleranceWellRate = 5e-3;

% Reset just in case
model.FlowDiscretization = [];
model.FlowPropertyFunctions = [];
model.PVTPropertyFunctions = [];
% Best performance
model.AutoDiffBackend.useMex = true;
model.AutoDiffBackend.rowMajor = true;

model = model.validateModel();
% Use the alternative more rigorous crossflow definition for component
% fluxes
xflow = WellComponentTotalVolumeBalanceCrossflow(model);
xflow.onlyLocalDerivatives = true;
model.FacilityModel.FacilityFlowDiscretization.PhaseFlux.allowCrossFlow = false;
model.FacilityModel.FacilityFlowDiscretization.ComponentTotalFlux = xflow;
% Disable flag for interpolation - better behavior
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
    'Name', 'GenericBlackOil_EGRID', 'nonlinearsolver', nls);

%% MRST grid     
simulatePackedProblem(problem, 'continueOnError', false);

%% Get results
[ws, states, reports] = getPackedSimulatorOutput(problem);
nstep = numel(ws);
T = cumsum(schedule.step.val(1:nstep));
%% Plot well sols
% We compare with different simulators. The biggest difference is that some
% of the simulations have hysteresis enabled. For the case without
% hystersis, MRST has excellent match.
ws_opm = output.opm.wellSols(2:end);
T_opm = cumsum(schedule.step.val);
fname = 'Flow legacy (hysteresis)';
wellSols = {ws, ws_opm};
time = {T; T_opm};
names = {'MRST (No hysteresis)', fname};

try
    % Flow output - may not be present
    pp = fullfile(mrstPath('opm-tests'), 'norne', 'NORNE_ATW2013');
    states_newopm = convertRestartToStates(pp, model.G);
    [ws_opm, T_opm_new] = convertSummaryToWellSols(pp);
    fname = 'Flow (Hysteresis)';
    
    wellSols{end+1} = ws_opm;
    time{end+1} = T_opm_new;
    names{end+1} = fname;
catch
    % Output not found
end

ecl = fullfile(mrstPath('opm-tests'), 'norne', 'ECL.2014.2', 'NORNE_ATW2013');
[ws_ecl, T_ecl] = convertSummaryToWellSols(ecl, 'metric');
wellSols{end+1} = ws_ecl;
time{end+1} = T_ecl;
names{end+1} = 'Eclipse (Hysteresis)';
fn = fullfile(mrstDataDirectory(), 'smry_nohyst.mat');

if ~exist(fn, 'file')
    url = 'https://www.sintef.no/contentassets/124f261f170947a6bc51dd76aea66129/smry_nohyst.mat';
    dfcn = mrstWebSave();
    dfcn(fn, url);
end
smry = load(fn);
smry = smry.smry;
[ws_nohyst, T_nohyst] = convertSummaryToWellSols(smry, 'metric');
wellSols{end+1} = ws_nohyst;
time{end+1} = T_nohyst;
names{end+1} = 'Flow (No hystersis)';

plotWellSols(wellSols, time, 'datasetnames', names, 'linestyles', {'-o', '--'})
%% Plot states
figure;
plotToolbar(G_viz, output.opm.states);
title('OPM');
figure;
plotToolbar(G_viz, states);
title('MRST');

%% Plot oil and gas rates
dt = schedule.step.val;

figure('position', 1e3*[ 1.0000    1.0977    0.8563    0.2403]);
hold on
for i = 1:2
    if i == 1
        wd = ws;
        t = T;
    else
        wd = ws_nohyst;
        t = T_nohyst;
    end
    dti = diff([0; t]);
    if isfield(wd{1}, 'status')
        status = getWellOutput(wd, 'status');
    else
        status = true;
    end
    getRates = @(x) max(-getWellOutput(wd, x).*status, 0);
    
    qo = getRates('qOs');
    qg = getRates('qGs');
    qw = getRates('qWs');

    qot = sum(qo, 2)/stb;
    qwt = sum(qw, 2)/stb;
    qgt = sum(qg, 2)/(1e3*ft^3);
    
    if i == 1
        plotarg = {'Linewidth', 2};
        style = '--';
    else
        plotarg = {};
        style = '-';
    end
    t = t/day;
    % Uncomment next line to get the cumulati
    % fn = @(x) cumsum(x.*dti, 1);
    fn = @(x) x;
    e = 0.5;
    plot(t, fn(qot), style, 'color', [1, e, e], plotarg{:});
    plot(t, fn(qwt), style, 'color', [e, e, 1], plotarg{:});
    plot(t, fn(qgt), style, 'color', [e, 1, e], plotarg{:});
    legend('Oil (stb)', 'Water (stb)', 'Gas (mscf)', 'location', 'northwest')
    axis tight
    xlabel('Time (days)');
end
%% Plot model with wells
k = rock.perm(:, 1);
figure;
plotCellData(G_viz, log10(k))
plotWell(G_viz, schedule.control(1).W, 'color', 'r', 'color2', 'b', 'FontSize', 8);
view(70, 30);
mrstColorbar(gca, k, 'southoutside', true);
axis tight

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>

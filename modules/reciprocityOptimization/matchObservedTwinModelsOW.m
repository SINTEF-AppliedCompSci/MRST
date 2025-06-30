function obj = matchObservedTwinModelsOW(model_1, states_1, schedule_1, model_2, states_2, schedule_2, observed, prior_dist, varargin)
% Compute mismatch between twin reservoir models with improved organization and scaling

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.
This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
%}

% Enhanced options with better default scaling
opt = struct('WaterRateWeight',     [], ...
    'OilRateWeight',       [], ...
    'BHPWeight',           [], ...
    'ComputePartials',     false, ...
    'tStep',              [], ...
    'state1',             [], ...
    'state2',             [], ...
    'from_states',        false, ...
    'matchOnlyProducers', false, ...
    'mismatchSum',        true, ...
    'accumulateWells',    [], ...
    'accumulateTypes',    [], ...
    'energyNormWeight',   0.0, ...  % Default regularization
    'H1NormWeight',       0.0, ...
    'L2NormWeight',       0.0, ...
    'interWellWeight',    0.00, ...
    'scalingMethod',      'dynamic'); % 'dynamic' or 'fixed'

opt = merge_options(opt, varargin{:});

% Time step handling
dts = schedule_1.step.val;
totTime = sum(dts);

if isempty(opt.tStep)
    tSteps = (1:numel(dts))';
else
    tSteps = opt.tStep;
    dts = dts(opt.tStep);
end

obj = repmat({[]}, numel(tSteps), 1);

for step = 1:numel(tSteps)
    current_step = tSteps(step);
    dt = dts(step);

    % Get all data in a more organized way
    [modelData1, modelData2, obsData] = getAllModelData(...
        model_1, states_1, schedule_1, ...
        model_2, states_2, schedule_2, ...
        observed, current_step, opt);

    % Calculate weights with improved scaling method
    %  [ww, wo, wp] = calculateWeights(...
    %      modelData1.qWs, modelData1.qOs, modelData1.bhp, ...
    %      modelData2.qWs, modelData2.qOs, modelData2.bhp, ...
    %         obsData.qWs, obsData.qOs, obsData.bhp, opt);
    [ww, wo, wp] = getWeights(obsData.qWs, obsData.qOs, obsData.bhp, opt);
    % Calculate all norms in separate organized functions
    norms = calculateNorms(...
        model_1, modelData1, model_2, modelData2, ...
        schedule_1, dt, totTime, opt);

    % Set normalization options
opt.normalizeNorms = true;
opt.useMobilityWeighting = true;

% Calculate norms
norms = calculateNormss(...
        model_1, modelData1, model_2, modelData2, ...
        schedule_1, dt, totTime, opt);

% Get weights considering mobility
% weights = getNormWeights(ref_model, {state1, state2}, opt);

    % Compute mismatch terms with proper scaling
    if opt.mismatchSum
        obj{step} = computeSummedMismatch(...
            modelData1, modelData2, obsData, ...
            ww, wo, wp, norms, dt, totTime, opt);
        if ~isempty(prior_dist)
            obj{step} = obj{step} + 0.0e-10.*sum(prior_dist{1}.^2);
        end
           
    else
        obj{step} = computeComponentMismatch(...
            modelData1, modelData2, obsData, ...
            ww, wo, wp, norms, dt, totTime, opt);
           
        if ~isempty(prior_dist)
            obj{step} = obj{step} + 0.0e-10.*sum(prior_dist{1}.^2);
        end
    end
end
end

%% Helper Functions

function [model1, model2, obs] = getAllModelData(model_1, states_1, schedule_1, ...
    model_2, states_2, schedule_2, ...
    observed, current_step, opt)
% Organized data extraction for both models and observed data

% Observed data
sol_obs = observed{current_step};
nw = numel(sol_obs.wellSol);
obs.status = vertcat(sol_obs.wellSol.status);
obs.qWs = vertcatIfPresent(sol_obs.wellSol, 'qWs', nw);
obs.qOs = vertcatIfPresent(sol_obs.wellSol, 'qOs', nw);
obs.bhp = vertcatIfPresent(sol_obs.wellSol, 'bhp', nw);
obs.sign = vertcat(sol_obs.wellSol.sign);

% Model 1 data
if opt.ComputePartials && opt.from_states
    state_1 = model_1.getStateAD(states_1{current_step}, true);
elseif opt.ComputePartials
    state_1 = opt.state1;
else
    state_1 = states_1{current_step};
end

model1 = extractModelData(model_1, state_1, true);

% Model 2 data
if opt.ComputePartials && opt.from_states
    state_2 = model_2.getStateAD(states_2{current_step}, true);
elseif opt.ComputePartials
    state_2 = opt.state2;
else
    state_2 = states_2{current_step};
end

model2 = extractModelData(model_2, state_2, true);

% Handle inactive wells consistently
[model1, model2, obs] = handleInactiveWells(model1, model2, obs);
end

function modelData = extractModelData(model, state, isAD)
% Extract all relevant data from a model state
nw = numel(state.wellSol);

modelData.status = vertcat(state.wellSol.status);

if isAD && isa(state.pressure, 'ADI')
    modelData.pressure = model.getProp(state, 'pressure');
    modelData.saturation = model.getProp(state, 's');
    modelData.qWs = model.FacilityModel.getProp(state, 'qWs');
    modelData.qOs = model.FacilityModel.getProp(state, 'qOs');
    modelData.bhp = model.FacilityModel.getProp(state, 'bhp');
    modelData.wellSol = state.wellSol;  % Store the entire wellSol structure
    modelData.totalFlux = state.FluxDisc.ComponentTotalFlux;
    modelData.totalMass = state.FlowProps.ComponentTotalMass;
    modelData.Mobility = state.FlowProps.Mobility;
    modelData.FaceComponentMobility = state.FluxDisc.FaceComponentMobility;

else
    modelData.pressure = state.pressure;
    modelData.saturation = state.s;
    modelData.totalFlux = state.FluxDisc.ComponentTotalFlux;
    modelData.totalMass = state.FlowProps.ComponentTotalMass;
    modelData.Mobility = state.FlowProps.Mobility;
    modelData.FaceComponentMobility = state.FluxDisc.FaceComponentMobility;
    modelData.qWs = vertcatIfPresent(state.wellSol, 'qWs', nw);
    modelData.qOs = vertcatIfPresent(state.wellSol, 'qOs', nw);
    modelData.bhp = vertcatIfPresent(state.wellSol, 'bhp', nw);
    modelData.wellSol = state.wellSol;  % Store the entire wellSol structure
end

% Ensure pressure is cell array for consistency
if ~iscell(modelData.pressure)
    modelData.pressure = {modelData.pressure};
    modelData.saturation = {modelData.saturation};
end
end

function [model1, model2, obs] = handleInactiveWells(model1, model2, obs)
% Handle inactive wells consistently across models and observations
if ~all(model1.status) || ~all(model2.status) || ~all(obs.status)
    [model1.bhp, obs.bhp] = expandToFull(model1.bhp, obs.bhp, model1.status, obs.status, true);
    [model1.qWs, obs.qWs] = expandToFull(model1.qWs, obs.qWs, model1.status, obs.status, false);
    [model1.qOs, obs.qOs] = expandToFull(model1.qOs, obs.qOs, model1.status, obs.status, false);

    [model2.bhp, obs.bhp] = expandToFull(model2.bhp, obs.bhp, model2.status, obs.status, true);
    [model2.qWs, obs.qWs] = expandToFull(model2.qWs, obs.qWs, model2.status, obs.status, false);
    [model2.qOs, obs.qOs] = expandToFull(model2.qOs, obs.qOs, model2.status, obs.status, false);
end
end

function [ww, wo, wp] = calculateWeights(qWs1, qOs1, bhp1, qWs2, qOs2, bhp2, qWs_obs, qOs_obs, bhp_obs, opt)
% Improved weight calculation with multiple methods
% Safe extraction of values (works for both ADI and numeric types)
qWs1 = getValue(qWs1);
qOs1 = getValue(qOs1);
bhp1 = getValue(bhp1);
qWs2 = getValue(qWs2);
qOs2 = getValue(qOs2);
bhp2 = getValue(bhp2);
switch opt.scalingMethod
    case 'dynamic'
        % Dynamic scaling based on observed and simulated values
        ww = 1 ./ max(mean(abs([qWs1, qWs2, qWs_obs])), eps);
        wo = 1 ./ max(mean(abs([qOs1, qOs2, qOs_obs])), eps);
        wp = 1 ./ max(mean(abs([bhp1, bhp2, bhp_obs])), eps);

    otherwise % 'fixed' or manual weights
        if isempty(opt.WaterRateWeight)
            ww = 1 ./ max(mean(abs(qWs_obs)), eps);
        else
            ww = opt.WaterRateWeight;
        end

        if isempty(opt.OilRateWeight)
            wo = 1 ./ max(mean(abs(qOs_obs)), eps);
        else
            wo = opt.OilRateWeight;
        end

        if isempty(opt.BHPWeight)
            wp = 1 ./ max(std(bhp_obs), eps);
        else
            wp = opt.BHPWeight;
        end
end

% Normalize weights to prevent any one term from dominating
total_weight = ww + wo + wp;
ww = ww / total_weight;
wo = wo / total_weight;
wp = wp / total_weight;
end

function val = getValue(x)
% Safely extracts value from ADI or numeric variables
if isa(x, 'ADI')
    val = value(x);
else
    val = x;
end
end

function  [ww, wo, wp] = getWeights(qWs, qOs, bhp, opt)
ww = opt.WaterRateWeight;
wo = opt.OilRateWeight;
wp = opt.BHPWeight;

rw = sum(abs(qWs)+abs(qOs));

if isempty(ww)
    % set to zero if all are zero
    if sum(abs(qWs))==0
        ww = 0;
    else
        ww = 1/rw;
    end
end

if isempty(wo)
    % set to zero if all are zero
    if sum(abs(qOs))==0
        wo = 0;
    else
        wo = 1/rw;
    end
end

if isempty(wp)
    % set to zero all are same
    dp = max(bhp)-min(bhp);
    if dp == 0
        wp = 0;
    else
        wp = 1/dp;
    end
end
end
function norms = calculateNorms(model_1, model1, model_2, model2, schedule, dt, totTime, opt)
% Calculate all norm terms with proper scaling

norms.energy = 0;
norms.H1 = 0;
norms.L2 = 0;
norms.interWell = 0;

% Energy Norm (Global Domain)
Div = model_1.operators.Div;
for i = 1:length(model1.totalMass)
    massDiff = (model1.totalMass{i} - model2.totalMass{i}) ./ dt;
    fluxDivDiff = Div(model1.totalFlux{i} - model2.totalFlux{i});
    norms.energy = norms.energy + (dt/totTime) * (sum(fluxDivDiff.^2) + sum(massDiff.^2));
end

% H1 Norm for Pressure
Grad = model_1.operators.Grad;
pressureDiff = model1.pressure{1} - model2.pressure{1};
gradPressure = Grad(model2.pressure{1});
norms.H1 = (dt/totTime) * sum(Grad(pressureDiff).^2) ./ max(sum(gradPressure.^2), eps);

% L2 Norm for Pressure at Wells
wellCells = reshape([schedule.control(1).W.cells], [], 1);
wellPressureDiff = (model1.pressure{1}(wellCells) - model2.pressure{1}(wellCells));
if strcmpi(opt.scalingMethod, 'dynamic')
    wellPressureDiff = wellPressureDiff ./ max(model2.pressure{1}(wellCells), eps);
end
norms.L2 = (dt/totTime) * sum(wellPressureDiff.^2);

% Inter-well Energy Norm
nw = numel(schedule.control(1).W);
for wew = 1:nw
    wc = schedule.control(1).W(wew).cells;
    wellPressureDiff = (model1.pressure{1}(wc) - model2.pressure{1}(wc));
    if strcmpi(opt.scalingMethod, 'dynamic')
        wellPressureDiff = wellPressureDiff ./ max(model2.pressure{1}(wc), eps);
    end

    % Access ComponentTotalFlux correctly from wellSol structure
    if isfield(model1.wellSol(wew), 'ComponentTotalFlux')
        flux1 = model1.wellSol(wew).ComponentTotalFlux;
        flux2 = model2.wellSol(wew).ComponentTotalFlux;

        if size(flux1, 2) > 0  % Check if flux data exists
            for i = 1:size(flux1, 2)
                norms.interWell = norms.interWell + ...
                    (dt/totTime) * sum(abs((flux1(:,i) - flux2(:,i)) .* wellPressureDiff));
            end
        end
    end
end
end

function totalMismatch = computeSummedMismatch(model1, model2, obs, ww, wo, wp, norms, dt, totTime, opt)
% Compute summed mismatch with all terms

matchCases = true(numel(obs.qWs), 1);
if opt.matchOnlyProducers
    matchCases = (obs.sign < 0);
end

time_scale = dt / (totTime * nnz(matchCases));

% Well terms
wellTerms = time_scale * sum(...
    (ww.*matchCases.*(model1.qWs - model2.qWs)).^2 + ...
    (wo.*matchCases.*(model1.qOs - model2.qOs)).^2 + ...
    (wp.*matchCases.*(model1.bhp - model2.bhp)).^2);

% Norm terms
normTerms = opt.energyNormWeight * norms.energy + ...
    opt.H1NormWeight * norms.H1 + ...
    opt.L2NormWeight * norms.energy + ...
    wp.*norms.well_pressure;

totalMismatch = wellTerms + normTerms;
end

function componentMismatch = computeComponentMismatch(model1, model2, obs, ww, wo, wp, norms, dt, totTime, opt)
% Compute individual mismatch components

matchCases = true(numel(obs.qWs), 1);
if opt.matchOnlyProducers
    matchCases = (obs.sign < 0);
end

fac = dt / (totTime * nnz(matchCases));

% Well components
wellComponents = {...
    fac * (ww.*matchCases.*(model1.qWs - obs.qWs)).^2, ...
    fac * (wo.*matchCases.*(model1.qOs - obs.qOs)).^2, ...
    fac * (wp.*matchCases.*(model1.bhp - obs.bhp)).^2, ...
    fac * (ww.*matchCases.*(obs.qWs - model2.qWs)).^2, ...
    fac * (wo.*matchCases.*(obs.qOs - model2.qOs)).^2, ...
    fac * (wp.*matchCases.*(obs.bhp - model2.bhp)).^2};


% Norm components
normComponents = {...
    opt.energyNormWeight * norms.energy, ...
    opt.H1NormWeight * norms.H1, ...
    wp.* norms.well_pressure, ...
    wp.* norms.well_pressure};

% Handle accumulation options
if isempty(opt.accumulateTypes)
    tmp = [wellComponents, normComponents];
else
    pt = opt.accumulateTypes;
    tmp = num2cell(zeros(1, max(pt)));
    for k = 1:numel(wellComponents)
        if pt(k) > 0
            tmp{pt(k)} = tmp{pt(k)} + wellComponents{k};
        end
    end
end

if ~isempty(opt.accumulateWells)
    pw = opt.accumulateWells;
    M = sparse(pw(pw>0), find(pw), 1);
    tmp = cellfun(@(x) M*x, tmp, 'UniformOutput', false);
end

componentMismatch = vertcat(tmp{:});
end

function v = vertcatIfPresent(sol, fn, nw)
if isfield(sol, fn)
    v = vertcat(sol.(fn));
    assert(numel(v) == nw);
    v = v(vertcat(sol.status));
else
    v = zeros(nnz(sol.status), 1);
end
end

function [v, v_obs] = expandToFull(v, v_obs, status, status_obs, setToZero)
tmp = zeros(size(status));
if isa(v, 'ADI')
    tmp = double2ADI(tmp, v);
end
tmp(status) = v;
v = tmp;

tmp = zeros(size(status));
tmp(status_obs) = v_obs;
v_obs = tmp;

if setToZero
    ix = status ~= status_obs;
    v(ix) = 0;
    v_obs(ix) = 0;
end
end

function norms = calculateNormss(model_1, model1, model_2, model2, schedule, dt, totTime, opt)
% Calculate comprehensive norms for reservoir optimization with mobility weighting
%
% Includes:
% 1. H1 norm of pressure differences (gradient term)
% 2. Energy norm from mass conservation (∂m/∂t + ∇·flux)
% 3. Well pressure mismatch terms
% 4. Flux continuity norm
% 5. Mobility-weighted versions of all terms

    % Initialize norm structure
    norms = struct('H1', 0, ...             % H1 pressure norm
                  'energy', 0, ...          % Mass conservation energy
                  'well_pressure', 0, ...   % Well pressure differences
                  'flux_continuity', 0, ... % Flux differences
                  'H1_mob', 0, ...         % Mobility-weighted H1
                  'energy_mob', 0);        % Mobility-weighted energy

    % Extract model operators and properties
    rock = model_1.rock;
    pv = model_1.operators.pv;
    T = model_1.operators.T;

    G = model_1.G;
    Div = model_1.operators.Div;
    Grad = model_1.operators.Grad;
    
    % Get pressure fields and differences
    p1 = model1.pressure{1};
    p2 = model2.pressure{1};
    dp = p1 - p2;
    faceMob1 =   model1.FaceComponentMobility;
    faceMob2 =   model2.FaceComponentMobility;
    %% 1. H1 Norm (Pressure Gradient Norm)
    % ----------------------------------
    norms.H1 = 0;
    for i = 1:sum(size(faceMob2))
        if ~isempty(faceMob1{i})
            grad_dp = faceMob1{i}.*T.*Grad(model1.pressure{1})-faceMob2{i}.*T.*Grad(model2.pressure{1});            
            norms.H1 = norms.H1 + (dt/totTime) * sum((grad_dp).^2);
        end
    end
    %% 2. Energy Norm (Mass Conservation)
    % ----------------------------------
    for i = 1:length(model1.totalMass)
        % Time derivative of mass
        dm_dt = (model1.totalMass{i} - model2.totalMass{i}) / dt;
        
        % Flux divergence difference
        div_flux_diff = Div(model1.totalFlux{i} - model2.totalFlux{i});
        
        % Combined energy term
        energy_term = dm_dt + div_flux_diff;
        norms.energy = norms.energy + (dt/totTime) * sum(pv .* energy_term.^2);
    end
    
    %% 3. Well Pressure Mismatch
    % --------------------------
    well_cells = vertcat(schedule.control(1).W.cells);
    if ~isempty(well_cells)
        dp_well = dp(well_cells);
        norms.well_pressure = (dt/totTime) * sum(pv(well_cells) .* dp_well.^2);
    end
  
    
   
    
   
end

function weights = getNormWeights(model, states, opt)
% Returns physically scaled weights using values from opt.weighting

    % Extract characteristic scales from weighting options
    if isfield(opt, 'weighting')
        % Get inverse of weights to estimate characteristic scales
        if isfield(opt.weighting, 'WaterRateWeight') && opt.weighting.WaterRateWeight > 0
            charQ = 1/sqrt(opt.weighting.WaterRateWeight);
        else
            charQ = 100*stb/day; % Default if not specified
        end
        
        if isfield(opt.weighting, 'BHPWeight') && opt.weighting.BHPWeight > 0
            charP = 1/sqrt(opt.weighting.BHPWeight);
        else
            charP = mean(states{1}.pressure); % Fallback to average pressure
        end
        
        % Estimate energy scale from pressure and flow rate
        charE = charP * charQ; % barsa * stb/day
    else
        % Fallback defaults
        charP = mean(states{1}.pressure);
        charQ = 100*stb/day;
        charE = charP * charQ;
    end

    weights = struct();
    
    % Base weights derived from characteristic scales
    weights.H1 = 1/(charP^2 * mean(model.rock.perm)); % H1 pressure norm
    weights.energy = 1/charE^2;                      % Energy norm
    weights.well_pressure = 1/charP^2;               % Well pressure
    weights.flux = 1/charQ^2;                        % Flux continuity
    
    % Apply user overrides if specified
    if isfield(opt, 'H1Weight') && opt.H1Weight > 0
        weights.H1 = opt.H1Weight;
    end
    if isfield(opt, 'energyWeight') && opt.energyWeight > 0
        weights.energy = opt.energyWeight;
    end
    
    % Mobility-weighted versions (if mobility data exists)
    if isfield(states{1}.FlowProps, 'Mobility')
        avg_mob = mean(states{1}.FlowProps.Mobility(:));
        weights.H1_mob = weights.H1 * avg_mob;
        weights.energy_mob = weights.energy * avg_mob;
    end
    
    % Normalize weights to [0,1] range while preserving ratios
    all_weights = struct2array(weights);
    total_weight = sum(all_weights(all_weights > 0)); % Only sum positive weights
    
    if total_weight > 0
        fn = fieldnames(weights);
        for i = 1:numel(fn)
            if weights.(fn{i}) > 0
                weights.(fn{i}) = weights.(fn{i}) / total_weight;
            end
        end
    end
    
    % Ensure minimum weighting for critical terms
    min_weight = 1e-4;
    if weights.energy < min_weight
        weights.energy = min_weight;
    end
    if weights.H1 < min_weight
        weights.H1 = min_weight; 
    end
end
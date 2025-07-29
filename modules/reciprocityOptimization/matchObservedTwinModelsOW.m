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
        model_1, states_1, ...
        model_2, states_2, ...
        observed, current_step, opt);

    % Calculate weights with improved scaling method
    [ww, wo, wp] = getWeights(obsData.qWs, obsData.qOs, obsData.bhp, opt);
    % Calculate all norms in separate organized functions


    % Set normalization options
    opt.normalizeNorms = true;
    opt.useMobilityWeighting = true;

    % Calculate norms
    norms = calculateNorms(...
        model_1, modelData1, modelData2, ...
        dt, totTime, opt);

    % Get weights considering mobility

    % Compute mismatch terms with proper scaling
    if opt.mismatchSum
        obj{step} = computeSummedMismatch(...
            modelData1, modelData2, schedule_1, obsData,...
            ww, wo, wp, norms, dt, totTime, opt);
        if ~isempty(prior_dist)
            obj{step} = obj{step} + 0.0e-10.*sum(prior_dist{1}.^2);
        end

    else
        obj{step} = computeComponentMismatch(...
            modelData1, modelData2, schedule_1, obsData , ...
            ww, wo, wp, norms, dt, totTime, opt);

        if ~isempty(prior_dist)
            obj{step} = obj{step} + 0.0e-10.*sum(prior_dist{1}.^2);
        end
    end
end
end

%% Helper Functions

function [model1, model2, obs] = getAllModelData(model_1, states_1, ...
    model_2, states_2, ...
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


function totalMismatch = computeSummedMismatch(model1, model2, schedule, obs, ww, wo, wp, norms, dt, totTime, opt)
% Compute summed mismatch with all terms

matchCases = true(numel(obs.qWs), 1);
if opt.matchOnlyProducers
    matchCases = (obs.sign < 0);
end

%  Well Pressure Mismatch
well_cells = vertcat(schedule.control(1).W.cells);
if ~isempty(well_cells)
    % Get pressure fields and differences
    p1 = model1.pressure{1};
    p2 = model2.pressure{1};
    dp = p1 - p2;
    dp_well = dp(well_cells);
end

fac = dt / (totTime * nnz(matchCases));

% Well terms
wellTerms = fac * sum(...
    (ww.*matchCases.*(model1.qWs - obs.qWs)).^2 + ...
    (wo.*matchCases.*(model1.qOs - obs.qOs)).^2 + ...
    (wp.*matchCases.*(model1.bhp - obs.bhp)).^2 + ...convMap
    (ww.*matchCases.*(obs.qWs - model2.qWs)).^2 + ...
    (wo.*matchCases.*(obs.qOs - model2.qOs)).^2 + ...
    (wp.*matchCases.*(obs.bhp - model2.bhp)).^2 + ...
    wp.*10*wp.*(matchCases.*(dp_well)).^2);

% Norm terms
normTerms = norms.energy + ...
      norms.H1 + 10*wp.*fac.*wp.*sum(((p1-p2)).^2);

totalMismatch = wellTerms + normTerms;
end

function componentMismatch = computeComponentMismatch(model1, model2, schedule, obs, ww, wo, wp, norms, dt, totTime, opt)
% Compute individual mismatch components

matchCases = true(numel(obs.qWs), 1);
if opt.matchOnlyProducers
    matchCases = (obs.sign < 0);
end
    
%  Well Pressure Mismatch
well_cells = vertcat(schedule.control(1).W.cells);
if ~isempty(well_cells)
    % Get pressure fields and differences
    p1 = model1.pressure{1};
    p2 = model2.pressure{1};
    dp = p1 - p2;
    dp_well = dp(well_cells);
end
fac = dt / (totTime * nnz(matchCases));

% Well components
wellComponents = {...
    fac * (ww.*matchCases.*(model1.qWs - obs.qWs)).^2, ...
    fac * (wo.*matchCases.*(model1.qOs - obs.qOs)).^2, ...
    fac * (wp.*matchCases.*(model1.bhp - obs.bhp)).^2, ...
    fac * (ww.*matchCases.*(obs.qWs - model2.qWs)).^2, ...
    fac * (wo.*matchCases.*(obs.qOs - model2.qOs)).^2, ...
    fac * (wp.*matchCases.*(obs.bhp - model2.bhp)).^2, ...
    0.3*fac * wp.*wp.*(matchCases.*(dp_well)).^2};



% Norm components
normComponents = {...
    0.001.*norms.energy, ...
    0.001.*norms.H1, 0.3*wp.*fac * wp.*sum(((p1-p2)).^2)};

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

function norms = calculateNorms(model_1, model1, model2, dt, totTime, opt)
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
                  'energy', 0);          % Mass conservation energy

    % Extract model operators and properties
    T = model_1.operators.T;

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
            grad_dp = faceMob1{i}.*T.*Grad(p1)-faceMob2{i}.*T.*Grad(p2);            
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
        norms.energy = norms.energy + (dt/totTime) * sum(energy_term.^2);
    end
  
   
end

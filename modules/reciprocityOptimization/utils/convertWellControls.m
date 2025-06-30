function newSchedule = convertWellControls(schedule, states, model, varargin)
% Convert well control types in a schedule and optionally add boundary conditions
%
% SYNOPSIS:
%   newSchedule = convertWellControls(schedule, states, model)
%   newSchedule = convertWellControls(schedule, states, model, 'pn', pv, ...)
%
% PARAMETERS:
%   schedule - Original MRST schedule structure
%   states   - Cell array of state structures for each timestep
%   model    - MRST reservoir model (for phase detection)
%
% OPTIONAL PARAMETERS:
%   'bc'            - Boundary conditions (cell array matching schedule steps)
%   'ConversionMap' - Struct specifying control type conversions
%   'Limits'        - Struct specifying well limits for converted controls
%   'Sign'          - Struct specifying well signs for converted controls

    % Parse input arguments
    opt = struct('bc', {{}}, ...
                 'ConversionMap', struct('bhp', 'rate', 'rate', 'bhp'), ...
                 'useLimits', false, ...
                 'Sign', struct());
    opt = merge_options(opt, varargin{:});
    
    % Get active phases from model (logical array: [water, oil, gas])
    activePhases = model.getActivePhases();
    phaseNames = {'water', 'oil', 'gas'};
    activePhaseNames = phaseNames(activePhases);
    
    % Validate conversion map based on active phases
    validateConversionMap(opt.ConversionMap, activePhases);
    
    % Initialize schedule
    nSteps = numel(schedule.step.val);
    nWells = numel(schedule.control(1).W);
    newScheduleSteps = cell(nSteps, 1);
    
    % Process each timestep
    for step = 1:nSteps
        currentControl = schedule.control(schedule.step.control(step));
        W = currentControl.W;
        
        % Convert each well's control type
        for w = 1:nWells
            wellSol = states{step}.wellSol(w);
            currentType = W(w).type;
            
            % Get conversion type for this well's current control type
            if isfield(opt.ConversionMap, currentType)
                newType = opt.ConversionMap.(currentType);
            else
                % Default to no conversion if type not specified
                newType = currentType;
            end
            
            % Perform the conversion
            [W(w), rateSum] = convertWellControl(W(w), wellSol, newType, activePhases);
            
            % Apply limits if specified
            W(w) = applyWellLimits(W(w), newType, opt.useLimits, wellSol, rateSum, activePhases);
            
            % Set well sign
            W(w) = setWellSign(W(w), newType, opt.Sign);
        end
        
        % Create schedule step with optional BCs
        if ~isempty(opt.bc) && numel(opt.bc) >= step && ~isempty(opt.bc{step})
            newScheduleSteps{step} = simpleSchedule(schedule.step.val(step), ...
                                                  'W', W, 'bc', opt.bc{step});
        else
            newScheduleSteps{step} = simpleSchedule(schedule.step.val(step), 'W', W);
        end
    end
    
    % Combine all steps into final schedule
    newSchedule = combineSchedules(newScheduleSteps{:}, 'makeConsistent', false);
end

%% Helper functions
function validateConversionMap(convMap, activePhases)
% Validate the conversion map structure based on active phases
    validFromTypes = {'bhp', 'rate'};
    
    % Determine valid target types based on active phases
    validToTypes = {'bhp'}; % Always valid
    if any(activePhases) % At least one phase active
        validToTypes{end+1} = 'rate';
    end
    if activePhases(2) % Oil active
        validToTypes{end+1} = 'orat';
    end
    if activePhases(1) % Water active
        validToTypes{end+1} = 'wrat';
    end
    if activePhases(3) % Gas active
        validToTypes{end+1} = 'grat';
    end
    
    fields = fieldnames(convMap);
    for i = 1:numel(fields)
        if ~ismember(fields{i}, validFromTypes)
            error('Invalid source control type: %s', fields{i});
        end
        if ~ismember(convMap.(fields{i}), validToTypes)
            error('Invalid target control type %s for active phases', ...
                  convMap.(fields{i}));
        end
    end
end

function [well, rateSum] = convertWellControl(well, wellSol, newType, activePhases)
% Convert a single well's control type with phase awareness
    phaseNames = {'water', 'oil', 'gas'};
    activePhaseNames = phaseNames(activePhases);
    
    rateSum = getRateSum(wellSol, activePhases);
    
    switch well.type
        case 'bhp'
            % Convert from BHP to rate control
            well.type = newType;
            switch newType
                case 'rate'
                    well.val = rateSum;
                    well.compi = getNormalizedCompi(wellSol, activePhases);
                case 'orat'
                    if ~activePhases(2)
                        error('Cannot convert to oil rate - oil phase not active');
                    end
                    well.val = getFieldOrDefault(wellSol, 'qOs', 0);
                    well.compi = setCompiForActivePhases(activePhases, 'oil');
                case 'wrat'
                    if ~activePhases(1)
                        error('Cannot convert to water rate - water phase not active');
                    end
                    well.val = getFieldOrDefault(wellSol, 'qWs', 0);
                    well.compi = setCompiForActivePhases(activePhases, 'water');
                case 'grat'
                    if ~activePhases(3)
                        error('Cannot convert to gas rate - gas phase not active');
                    end
                    well.val = getFieldOrDefault(wellSol, 'qGs', 0);
                    well.compi = setCompiForActivePhases(activePhases, 'gas');
            end
            
        case 'rate'
            % Convert from rate to BHP control
            if strcmp(newType, 'bhp')
                well.type = 'bhp';
                well.val = wellSol.bhp;
                well.compi = getNormalizedCompi(wellSol, activePhases);
            end
    end
    
    % Clear any existing limits
    well.lims = [];
end

function rateSum = getRateSum(wellSol, activePhases)
% Calculate total surface rate sum based on active phases
    rateSum = 0;
    
    if activePhases(1) % Water
        rateSum = rateSum + (getFieldOrDefault(wellSol, 'qWs', 0));
    end
    if activePhases(2) % Oil
        rateSum = rateSum + (getFieldOrDefault(wellSol, 'qOs', 0));
    end
    if activePhases(3) % Gas
        rateSum = rateSum + (getFieldOrDefault(wellSol, 'qGs', 0));
    end
end

function compi = getNormalizedCompi(wellSol, activePhases)
% Get normalized compi from wellSol rates based on active phases
    compi = zeros(1,sum(activePhases)); % [water, oil, gas]
    
    if activePhases(1)
        compi(1) = abs(getFieldOrDefault(wellSol, 'qWs', 0));
    end
    if activePhases(2)
        compi(2) = abs(getFieldOrDefault(wellSol, 'qOs', 0));
    end
    if activePhases(3)
        compi(3) = abs(getFieldOrDefault(wellSol, 'qGs', 0));
    end
    
    rateSum = sum(compi);
    
    if (rateSum~=0)
        compi = compi / rateSum;
    end
    
    % Return only compi for active phases
    compi = compi(activePhases);
end

function compi = setCompiForActivePhases(activePhases, targetPhase)
% Create compi vector for single-phase control
    compi = zeros(1,sum(activePhases ...
        )); % [water, oil, gas]
    
    switch targetPhase
        case 'water'
            compi(1) = 1;
        case 'oil'
            compi(2) = 1;
        case 'gas'
            compi(3) = 1;
    end
    
    % Return only compi for active phases
    compi = compi(activePhases);
end

function well = applyWellLimits(well, newType, useLimits, wellSol, rateSum, activePhases)
% Apply limits to converted well with phase awareness
    if ~useLimits
        return;
    end
    
    well.lims = struct();
    
    switch newType
        case 'bhp'
            if isfield(wellSol, 'qWs')
                well.lims.qWs = wellSol.qWs;
            end
               
            if isfield(wellSol, 'qOs')
                well.lims.qOs = wellSol.qOs;
            end
               
            if isfield(wellSol, 'qGs')
                well.lims.qGs = wellSol.qGs;
            end
            
        case {'rate', 'orat', 'wrat', 'grat'}
            if isfield(wellSol, 'bhp')
                well.lims.bhp = wellSol.bhp;
            end
            
    end
end

function well = setWellSign(well, newType, signStruct)
% Set well sign based on conversion type
    
end

function value = getFieldOrDefault(structure, fieldName, defaultValue)
% Get field value or return default if missing/empty
    if isfield(structure, fieldName) && ~isempty(structure.(fieldName))
        value = structure.(fieldName);
    else
        value = defaultValue;
    end
end
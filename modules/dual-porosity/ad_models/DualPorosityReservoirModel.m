classdef DualPorosityReservoirModel < PhysicalModel
    %Base class for physical models
    %
    % SYNOPSIS:
    %   model = ReservoirModel(G, rock, fluid)
    %
    % DESCRIPTION:
    %   Extension of `PhysicalModel` class to accomodate reservoir-specific
    %   features such as fluid and rock as well as commonly used phases and
    %   variables.
    %
    % REQUIRED PARAMETERS:
    %   G     - Simulation grid.
    %
    %   rock  - Valid rock used for the model.
    %
    %   fluid - Fluid model used for the model.
    %
    %
    % OPTIONAL PARAMETERS:
    %   'property' - Set property to the specified value.
    %
    % RETURNS:
    %   Class instance.
    %
    % SEE ALSO:
    %   `ThreePhaseBlackOilModel`, `TwoPhaseOilWaterModel`, `PhysicalModel`

properties
    % Required arguments
    fluid % The fluid model. See `initSimpleADIFluid`, `initDeckADIFluid`
    rock % The rock structure (perm/poro/ntg). See `makeRock`.
    % Limits
    dpMaxRel % Maximum relative pressure change
    dpMaxAbs % Maximum absolute pressure change
    dsMaxAbs % Maximum absolute saturation change
    maximumPressure % Maximum pressure allowed in reservoir
    minimumPressure % Minimum pressure allowed in reservoir
    % Phases and components
    water % Indicator showing if the aqueous/water phase is present
    gas % Indicator showing if the vapor/gas phase is present
    oil % Indicator showing if the liquid/oil phase is present
    % Tolerances
    useCNVConvergence % Use volume-scaled tolerance scheme
    toleranceCNV; % CNV tolerance (similar to inf-norm over saturation error)
    toleranceMB; % MB tolerance values (sum of mass-balance error)
    % Input/output
    inputdata % Input data used to instantiate the model
    extraStateOutput % Write extra data to states. Depends on submodel type.
    extraWellSolOutput % Output extra data to wellSols: GOR, WOR, ...
    outputFluxes % Store integrated fluxes in state.
    % Coupling to forces and other models
    gravity % Vector for the gravitational force
    FacilityModel % Facility model used to represent wells
    FlowPropertyFunctions % Grouping for flow properties
    FlowDiscretization % Grouping for flux discretization
    PVTPropertyFunctions
    Components = {};
    
    %% Dual Porosity properties
    transfer_model_object
    rock_matrix
    fluid_matrix
    
end

methods
    % --------------------------------------------------------------------%
    function model = DualPorosityReservoirModel(G, varargin)
        model = model@PhysicalModel(G);
        

        if nargin == 1 || ischar(varargin{1})
            % We were given only grid + any keyword arguments
            doSetup = false;
        else
            assert(nargin >= 3)
            rock_all = varargin{1};
            fluid_all = varargin{2};
            % We are being called in format
            % ReservoirModel(G, rock, fluid, ...)
            model.rock  = rock_all{1};
            model.fluid = fluid_all{1};

            % Rest of arguments should be keyword/value pairs.
            varargin = varargin(3:end);
            % We have been provided the means, so we will execute setup
            % phase after parsing other inputs and defaults.
            doSetup = ~(isempty(G) || isempty(model.rock));
        end
        
        %% Dual porosity: rock and fluid fracture/matrix
        model.rock_matrix = rock_all{2};
        model.fluid_matrix = fluid_all{2};
        
        %% Dual porosity: adding 'dummy' transfer object
        model.transfer_model_object = TransferFunction();
        
        %%
        model.dpMaxRel = inf;
        model.dpMaxAbs = inf;

        model.minimumPressure = -inf;
        model.maximumPressure =  inf;

        model.dsMaxAbs = .2;

        model.nonlinearTolerance = 1e-6;
        model.inputdata = [];

        model.useCNVConvergence = false;
        model.toleranceCNV = 1e-3;
        model.toleranceMB = inf;

        model.extraStateOutput = false;
        model.extraWellSolOutput = true;
        model.outputFluxes = true;
        % Gravity defaults to the global variable
        model.gravity = gravity();
        [model, unparsed] = merge_options(model, varargin{:}); %#ok

        % Base class does not support any phases
        model.water = false;
        model.gas = false;
        model.oil = false;

        if doSetup
            %% dual porosity: rock fracture here
            model.operators = setupOperatorsTPFA(G, model.rock, 'deck', model.inputdata);
            
            %% Dual porosity: matrix pore volume
            operators_matrix = setupOperatorsTPFA(G, model.rock_matrix, 'deck', model.inputdata);
            model.operators.pv_matrix = operators_matrix.pv;
        end
    end

    % --------------------------------------------------------------------%
    function state = validateState(model, state)
        % Validate initial state.
        %
        % SEE ALSO:
        %   :meth:`ad_core.models.PhysicalModel.validateState`
        state = validateState@PhysicalModel(model, state);
        active = model.getActivePhases();
        nPh = nnz(active);
        nc = model.G.cells.num;
        model.checkProperty(state, 'Pressure', [nc, 1], [1, 2]);
        if nPh > 1
            model.checkProperty(state, 'Saturation', [nc, nPh], [1, 2]);
        end
        if ~isempty(model.FacilityModel)
            state = model.FacilityModel.validateState(state);
        end
        if ~isfield(state, 'sMax') && isfield(state, 's')
            state.sMax = state.s;
        end
    end

    function vars = getSaturationVarNames(model)
        vars = {'sw', 'so', 'sg'};
        ph = model.getActivePhases();
        vars = vars(ph);
    end
    
    function vars = getSaturationVarNamesMatrix(model)
        vars = {'swm', 'som', 'sgm'};
        ph = model.getActivePhases();
        vars = vars(ph);
    end

    % --------------------------------------------------------------------%
    function dt = getMaximumTimestep(model, state, state0, dt, drivingForces)
        % Define the maximum allowable time-step based on physics or
        % discretization choice
        dt = getMaximumTimestep@PhysicalModel(model, state, state0, dt, drivingForces);
        if ~isempty(model.FlowDiscretization)
            dt = model.FlowDiscretization.getMaximumTimestep(model, state, state0, dt, drivingForces);
        end
    end

    function [model, state] = prepareReportstep(model, state, state0, dt, drivingForces)
        if ~isempty(drivingForces.W)
            assert(~isempty(model.FacilityModel), ...
            'FacilityModel not set up. Call model.validateModel before calling equations with wells.');
            [model.FacilityModel, state] = model.FacilityModel.prepareReportstep(state, state0, dt, drivingForces);
        end
    end

    function [model, state] = prepareTimestep(model, state, state0, dt, drivingForces)
        if ~isempty(drivingForces.W)
            [model.FacilityModel, state] = model.FacilityModel.prepareTimestep(state, state0, dt, drivingForces);
        end
        if ~isempty(model.FlowDiscretization)
            [model.FlowDiscretization, state] = model.FlowDiscretization.prepareTimestep(model, state, state0, dt, drivingForces);
        end
    end

    % --------------------------------------------------------------------%
    function model = validateModel(model, varargin)
        % Validate model.
        %
        % SEE ALSO:
        %   :meth:`ad_core.models.PhysicalModel.validateModel`
        model = validateModel@PhysicalModel(model, varargin{:});

        if isempty(model.FacilityModel)
            model.FacilityModel = FacilityModel(model); %#ok
        else
            model.FacilityModel.ReservoirModel = model;
        end

        assert(~isempty(model.operators),...
            'Operators must be set up before simulation. See model.setupOperators for more details.');

        model.FacilityModel = model.FacilityModel.validateModel(varargin{:});

        if isempty(model.FlowPropertyFunctions)
            model.FlowPropertyFunctions = FlowPropertyFunctions(model); %#ok
        end
        if isempty(model.FlowDiscretization)
            model.FlowDiscretization = FlowDiscretization(model); %#ok
        end
        if isempty(model.PVTPropertyFunctions)
            model.PVTPropertyFunctions = PVTPropertyFunctions(model); %#ok
        end
    end

    % --------------------------------------------------------------------%
    function [model, state] = updateForChangedControls(model, state, forces)
        % Called whenever controls change.
        %
        % NOTE:
        %   The addition this class makes is also updating the well
        %   solution and the underlying `FacilityModel` class instance.
        %
        % SEE ALSO:
        %   :meth:`ad_core.models.PhysicalModel.updateForChangedControls`

        % Called whenever controls change. Since this model can be used
        % with wells, we call the facility model's setup routine.
        model.FacilityModel = model.FacilityModel.setupWells(forces.W);
        state.wellSol = initWellSolAD(forces.W, model, state);
        [model, state] = updateForChangedControls@PhysicalModel(model, state, forces);
    end

    % --------------------------------------------------------------------%
    function [state, report] = updateState(model, state, problem, dx, drivingForces)
        % Generic update function for reservoir models containing wells.
        %
        % SEE ALSO:
        %   :meth:`ad_core.models.PhysicalModel.updateState`

        % Split variables into three categories: Regular/rest variables, saturation
        % variables (which sum to 1 after updates) and well variables (which live
        % in wellSol and are in general more messy to work with).
        [restVars, satVars, satMatrixVars, wellVars] = model.splitPrimaryVariables(problem.primaryVariables);

        % Update the wells
        if isfield(state, 'wellSol')
            state.wellSol = model.FacilityModel.updateWellSol(state.wellSol, problem, dx, drivingForces, wellVars);
        end

        % Update saturations in one go
        state  = model.updateSaturations(state, dx, problem, satVars);
        
        % Update saturations in one go
        state  = model.updateSaturationsMatrix(state, dx, problem, satMatrixVars);

        if ~isempty(restVars)
            % Handle pressure seperately
            state = model.updateStateFromIncrement(state, dx, problem, 'pressure', model.dpMaxRel, model.dpMaxAbs);
            state = model.capProperty(state, 'pressure', model.minimumPressure, model.maximumPressure);
            restVars = model.stripVars(restVars, 'pressure');

            % Update remaining variables (tracers, temperature etc)
            for i = 1:numel(restVars)
                 state = model.updateStateFromIncrement(state, dx, problem, restVars{i});
            end
        end

        report = [];
    end    
    
    % --------------------------------------------------------------------%
%     function [state, report] = updateState(model, state, problem, dx, drivingForces)
%         % Generic update function for reservoir models containing wells.
%         %
%         % SEE ALSO:
%         %   :meth:`ad_core.models.PhysicalModel.updateState`
% 
%         % Split variables into three categories: Regular/rest variables, saturation
%         % variables (which sum to 1 after updates) and well variables (which live
%         % in wellSol and are in general more messy to work with).
%         [restVars, satVars, wellVars] = model.splitPrimaryVariables(problem.primaryVariables);
% 
%         % Update the wells
%         if isfield(state, 'wellSol')
%             state.wellSol = model.FacilityModel.updateWellSol(state.wellSol, problem, dx, drivingForces, wellVars);
%         end
% 
%         % Update saturations in one go
%         state  = model.updateSaturations(state, dx, problem, satVars);
%         
%         % Update saturations in one go
% %         state  = model.updateSaturationsMatrix(state, dx, problem, satMatrixVars);
% 
%         if ~isempty(restVars)
%             % Handle pressure seperately
%             state = model.updateStateFromIncrement(state, dx, problem, 'pressure', model.dpMaxRel, model.dpMaxAbs);
%             state = model.capProperty(state, 'pressure', model.minimumPressure, model.maximumPressure);
%             restVars = model.stripVars(restVars, 'pressure');
% 
%             % Update remaining variables (tracers, temperature etc)
%             for i = 1:numel(restVars)
%                  state = model.updateStateFromIncrement(state, dx, problem, restVars{i});
%             end
%         end
% 
%         report = [];
%     end
    % --------------------------------------------------------------------%
    function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
        % Generic update function for reservoir models containing wells.
        %
        % SEE ALSO:
        %   :meth:`ad_core.models.PhysicalModel.updateAfterConvergence`
        if ~isempty(model.FacilityModel)
            [state, f_report] = model.FacilityModel.updateAfterConvergence(state0, state, dt, drivingForces);
        else
            f_report = [];
        end
        [state, report] = updateAfterConvergence@PhysicalModel(model, state0, state, dt, drivingForces);
        report.FacilityReport = f_report;
        if isfield(state, 'sMax')
            if ~all(size(state.sMax) == size(state.s))
                state.sMax = state.s;
            else
                state.sMax = max(state.sMax, state.s);
            end
        end
        if isfield(state, 'FacilityState')
            state = rmfield(state, 'FacilityState');
        end
    end

    % --------------------------------------------------------------------%
    function model = setupOperators(model, G, rock, varargin)
        % Set up default discretization operators and other static props
        %
        % SYNOPSIS:
        %   model = model.setupOperators(G, rock)
        %
        % DESCRIPTION:
        %   This function calls the default set of discrete operators as
        %   implemented by `setupOperatorsTPFA`. The default operators use
        %   standard choices for reservoir simulation similar to what is
        %   found in commercial simulators: Single-point potential upwind
        %   and two-point flux approximation.
        %
        % PARAMETERS:
        %   model - Class instance.
        %   G     - The grid used for computing the operators. Must have
        %           geometry information added (typically from
        %           `computeGeometry`). Although this is a property on the
        %           model, we allow for passing of a different grid for the
        %           operator setup since this is useful in some workflows.
        %   rock  - Rock structure. See `makeRock`.
        %
        % RETURNS:
        %   model - Model with updated `operators` property.
        %
        % NOTE:
        %   This function is called automatically during class
        %   construction.
        model.operators = setupOperatorsTPFA(G, rock, varargin{:});
    end

    % --------------------------------------------------------------------%

    function [values, tolerances, names] = getConvergenceValues(model, problem, varargin)
        % Check convergence criterion. Relies on `FacilityModel` to check
        % wells.
        %
        % SEE ALSO:
        %   :meth:`ad_core.models.PhysicalModel.checkConvergence`
        if model.useCNVConvergence
            % Use convergence model similar to commercial simulator
            [v_cnv, names_cnv, tol_cnv, is_cnv] = getConvergenceValuesCNV(model, problem);
            [v_wells, tol_wells, names_wells, is_well] = ...
                model.FacilityModel.getFacilityConvergenceValues(problem);
            % Get the values for all equations, just in case there are some
            % values that are not either wells or standard 3ph conservation
            % equations
            rest = ~(is_cnv | is_well);

            values_all = norm(problem, inf);
            names_rest = problem.equationNames(rest);
            values_rest = values_all(rest);
            tol_rest = repmat(model.nonlinearTolerance, size(values_rest));

            values = [v_cnv, v_wells, values_rest];
            tolerances = [tol_cnv, tol_wells, tol_rest];
            names = horzcat(names_cnv, names_wells, names_rest);
        else
            % Use strict tolerances on the residual without any
            % fingerspitzengefuhlen by calling the parent class
            [values, tolerances, names] = getConvergenceValues@PhysicalModel(model, problem, varargin{:});
            [v_wells, tol_wells, names_wells, is_well] = ...
                model.FacilityModel.getFacilityConvergenceValues(problem);
            if any(is_well)
                tolerances(is_well) = tol_wells;
                values(is_well) = v_wells;
                names(is_well) = names_wells;
            end
        end
    end

    % --------------------------------------------------------------------%
    function [fn, index] = getVariableField(model, name, varargin)
        % Map variables to state field.
        %
        % SEE ALSO:
        %   :meth:`ad_core.models.PhysicalModel.getVariableField`
        %% Dual porosity: variable names
        switch(lower(name))
            case {'t', 'temperature'}
                fn = 'T';
                index = ':';
            case {'swmax', 'somax', 'sgmax'}
                fn = 'sMax';
                index = model.satVarIndex(name(1:2));
            case 'smax'
                fn = 'sMax';
                index = ':';
            case {'sw', 'water', 'swf'}
                index = model.satVarIndex('sw');
                fn = 's';
            case {'so', 'oil', 'sof'}
                index = model.satVarIndex('so');
                fn = 's';
            case {'sg', 'gas', 'sgf'}
                index = model.satVarIndex('sg');
                fn = 's';
            case {'s', 'sat', 'saturation', 'sf'}
                index = ':';
                fn = 's';
            case {'pressure', 'p', 'pf'}
                index = ':';
                fn = 'pressure';
            case {'swm', 'water_matrix'}
                index = model.satVarIndexMatrix('swm');
                fn = 'sm';
            case {'som', 'oil_matrix'}
                index = model.satVarIndexMatrix('som');
                fn = 'sm';
            case {'sgm', 'gas_matrix'}
                index = model.satVarIndexMatrix('sgm');
                fn = 'sm';
            case {'sm', 'satm', 'saturation_matrix'}
                index = ':';
                fn = 'sm';
            case {'pressure_matrix', 'pm'}
                index = ':';
                fn = 'pressure_matrix';
            case 'wellsol'
                % Use colon to get all variables, since the wellsol may
                % be empty
                index = ':';
                fn = 'wellSol';
            otherwise
                % This will throw an error for us
                [fn, index] = getVariableField@PhysicalModel(model, name, varargin{:});
        end
    end

    function containers = getStateFunctionGroupings(model)
        containers = getStateFunctionGroupings@PhysicalModel(model);
        extra = {model.FlowPropertyFunctions, model.FlowDiscretization, model.PVTPropertyFunctions};
        if ~isempty(model.FacilityModel)
            fm_props = model.FacilityModel.getStateFunctionGroupings();
            extra = [extra, fm_props];
        end
        extra = extra(~cellfun(@isempty, extra));
        containers = [containers, extra];
    end
    % --------------------------------------------------------------------%
    function names = getComponentNames(model) %#ok
        % Get the names of components for the model
        % SYNOPSIS:
        %   names = model.getComponentNames();
        %
        % PARAMETERS:
        %   model - Class instance.
        %
        % RETURNS:
        %   names - Cell array of component names in no particular order.
        %
        names = {};
    end

    % --------------------------------------------------------------------%
    function forces = getValidDrivingForces(model)
        % Get valid forces. This class adds support for wells, bc and src.
        %
        % SEE ALSO:
        %   :meth:`ad_core.models.PhysicalModel.getValidDrivingForces`

        forces = getValidDrivingForces@PhysicalModel(model);
        % Support for wells
        forces.W   = [];
        % Support for boundary conditions
        forces.bc  = [];
        % Support for direct source terms
        forces.src = [];
    end

%     function [restVars, satVars, wellVars] = splitPrimaryVariables(model, vars)
%         % Split cell array of primary variables into grouping
%         % SYNOPSIS:
%         %   [restVars, satVars, wellVars] = model.splitPrimaryVariables(vars)
%         %
%         % DESCRIPTION:
%         %   Split a set of primary variables into three groups:
%         %   Well variables, saturation variables and the rest. This is
%         %   useful because the saturation variables usually are updated
%         %   together, and the well variables are a special case.
%         %
%         % PARAMETERS:
%         %   model - Class instance.
%         %   vars  - Cell array with names of primary variables
%         %
%         % RETURNS:
%         %   restVars - Names of variables that are not saturations or
%         %              belong to the wells.
%         %   satVars  - Names of the saturation variables present in `vars`.
%         %   wellVars - Names of the well variables present in `vars`
%         wellvars = model.FacilityModel.getPrimaryVariableNames();
%         isSat   = cellfun(@(x) any(strcmpi(model.getSaturationVarNames, x)), vars);
%         isWells = cellfun(@(x) any(strcmpi(wellvars, x)), vars);
% 
%         wellVars = vars(isWells);
%         satVars  = vars(isSat);
% 
%         restVars = vars(~isSat & ~isWells);
%     end
    
    
    % --------------------------------------------------------------------%
    function [restVars, satVars, satMatrixVars, wellVars] = splitPrimaryVariables(model, vars)
        % Split cell array of primary variables into grouping
        % SYNOPSIS:
        %   [restVars, satVars, wellVars] = model.splitPrimaryVariables(vars)
        %
        % DESCRIPTION:
        %   Split a set of primary variables into three groups:
        %   Well variables, saturation variables and the rest. This is
        %   useful because the saturation variables usually are updated
        %   together, and the well variables are a special case.
        %
        % PARAMETERS:
        %   model - Class instance.
        %   vars  - Cell array with names of primary variables
        %
        % RETURNS:
        %   restVars - Names of variables that are not saturations or
        %              belong to the wells.
        %   satVars  - Names of the saturation variables present in `vars`.
        %   wellVars - Names of the well variables present in `vars`
        wellvars = model.FacilityModel.getPrimaryVariableNames();
        isSat   = cellfun(@(x) any(strcmpi(model.getSaturationVarNames, x)), vars);
        isWells = cellfun(@(x) any(strcmpi(wellvars, x)), vars);
        isSatMatrix = cellfun(@(x) any(strcmpi(model.getSaturationVarNamesMatrix, x)), vars);

        wellVars = vars(isWells);
        satVars  = vars(isSat);
        satMatrixVars = vars(isSatMatrix);

        restVars = vars(~isSat & ~isWells & ~isSatMatrix);
    end

    % --------------------------------------------------------------------%
    function [isActive, phInd] = getActivePhases(model)
        % Get active flag for MRST's canonical phase ordering (WOG)
        %
        % SYNOPSIS:
        %   [act, indices] = model.getActivePhases();
        %
        % PARAMETERS:
        %   model    - Class instance
        %
        % RETURNS:
        %   isActive - Total number of known phases array with booleans
        %              indicating if that phase is present. MRST uses a
        %              ordering of water, oil and then gas.
        %
        %   phInd    - Indices of the phases present. For instance, if
        %              water and gas are the only ones present,
        %              `phInd = [1, 3]`

        isActive = [model.water, model.oil, model.gas];
        if nargout > 1
            phInd = find(isActive);
        end
    end

    % --------------------------------------------------------------------%
    function [phNames, longNames] = getPhaseNames(model)
        % Get short and long names of the present phases.
        %
        % SYNOPSIS:
        %   [phNames, longNames] = model.getPhaseNames();
        %
        % PARAMETERS:
        %   model    - Class instance
        %
        % RETURNS:
        %   phNames   - Cell array containing the short hanes ('W', 'O',
        %               G') of the phases present
        %
        %   longNames - Longer names ('water', 'oil', 'gas') of the phases
        %               present.
        tmp = 'WOG';
        active = model.getActivePhases();
        phNames = tmp(active);
        if nargout > 1
            tmp = {'water', 'oil', 'gas'};
            longNames = tmp(active);
        end
    end

    function phIndices = getPhaseIndices(model)
        % Get the active phases in canonical ordering
        w = model.water;
        o = model.oil;
        g = model.gas;
        phIndices = [w, w+o, w+o+g];
        phIndices(~model.getActivePhases) = -1;
    end

    % --------------------------------------------------------------------%
    function index = getPhaseIndex(model, phasename)
        % Query the index of a phase in the model's ordering
        %
        % SYNOPSIS:
        %   index = model.getPhaseNames();
        %
        % PARAMETERS:
        %   model     - Class instance
        %   phasename - The name of the phase to be found.
        % RETURNS:
        %   index   - Index of phase `phasename`
        active = model.getPhaseNames();
        index = find(active == phasename);
    end

    % --------------------------------------------------------------------%
    function n = getNumberOfComponents(model)
        n = numel(model.Components);
    end

    function n = getNumberOfPhases(model)
        n = sum(model.getActivePhases());
    end
    
    % --------------------------------------------------------------------%
    function state = updateSaturationsMatrix(model, state, dx, problem, satMatrixVars)
        % Update of phase-saturations
        %
        % SYNOPSIS:
        %   state = model.updateSaturations(state, dx, problem, satVars)
        %
        % DESCRIPTION:
        %   Update saturations (likely state.s) under the constraint that
        %   the sum of volume fractions is always equal to 1. This
        %   assumes that we have solved for n - 1 phases when n phases
        %   are present.
        %
        % PARAMETERS:
        %   model   - Class instance
        %   state   - State to be updated
        %   dx      - Cell array of increments, some of which correspond
        %             to saturations
        %   problem - `LinearizedProblemAD` class instance from which `dx`
        %             was obtained.
        %   satVars - Cell array with the names of the saturation
        %             variables.
        %
        % RETURNS:
        %   state - Updated state with saturations within physical
        %           constraints.
        %
        % SEE ALSO:
        %   `splitPrimaryVariables`

        if nargin < 5
            % Get the saturation names directly from the problem
            [~, satMatrixVars] = ...
                splitPrimaryVariables(model, problem.primaryVariables);
        end
        if isempty(satMatrixVars)
            % No saturations passed, nothing to do here.
            return
        end
        state_init = state;
        % Solution variables should be saturations directly, find the missing
        % link
        saturations = lower(model.getSaturationVarNamesMatrix);
        fillsat = setdiff(saturations, lower(satMatrixVars));

        n_fill = numel(fillsat);
        assert(n_fill == 1 || n_fill == 0)
        if n_fill == 1
            fillsat = fillsat{1};
            % Fill component is whichever saturation is assumed to fill up the rest of
            % the pores. This is done by setting that increment equal to the
            % negation of all others so that sum(s) == 0 at end of update
            solvedFor = ~strcmpi(saturations, fillsat);
        else
            % All saturations are primary variables. Sum of saturations is
            % assumed to be enforced from the equation setup
            solvedFor = true(numel(saturations), 1);
        end
        ds = zeros(model.G.cells.num, numel(saturations));

        tmp = 0;
        for i = 1:numel(saturations)
            if solvedFor(i)
                v = model.getIncrement(dx, problem, saturations{i});
                ds(:, i) = v;
                if n_fill > 0
                    % Saturations added for active variables must be subtracted
                    % from the last phase
                    tmp = tmp - v;
                end
            end
        end
        ds(:, ~solvedFor) = tmp;
        % We update all saturations simultanously, since this does not bias the
        % increment towards one phase in particular.
        state = model.updateStateFromIncrement(state, ds, problem, 'sm', inf, model.dsMaxAbs);
%         if isempty(model.FlowPropertyFunctions)
%             chopped = false(model.G.cells.num, 1);
%         else
%             kr = model.FlowPropertyFunctions.RelativePermeability;
%             [state, chopped] = kr.applyImmobileChop(model, state, state_init);
%         end
        if n_fill == 1
            % Ensure that values are within zero->one interval, and
            % re-normalize if any values were capped
            bad = any((state.sm > 1) | (state.sm < 0), 2);
            if any(bad)
                state.sm(bad, :) = min(state.sm(bad, :), 1);
                state.sm(bad, :) = max(state.sm(bad, :), 0);
                state.sm(bad, :) = bsxfun(@rdivide, state.sm(bad, :), sum(state.sm(bad, :), 2));
            end
        else
            % Ensure positive values, we assume that sum of saturations is
            % handled via constraint equation
            state.sm = max(state.sm, 0);
        end
    end

    % --------------------------------------------------------------------%
    function state = updateSaturations(model, state, dx, problem, satVars)
        % Update of phase-saturations
        %
        % SYNOPSIS:
        %   state = model.updateSaturations(state, dx, problem, satVars)
        %
        % DESCRIPTION:
        %   Update saturations (likely state.s) under the constraint that
        %   the sum of volume fractions is always equal to 1. This
        %   assumes that we have solved for n - 1 phases when n phases
        %   are present.
        %
        % PARAMETERS:
        %   model   - Class instance
        %   state   - State to be updated
        %   dx      - Cell array of increments, some of which correspond
        %             to saturations
        %   problem - `LinearizedProblemAD` class instance from which `dx`
        %             was obtained.
        %   satVars - Cell array with the names of the saturation
        %             variables.
        %
        % RETURNS:
        %   state - Updated state with saturations within physical
        %           constraints.
        %
        % SEE ALSO:
        %   `splitPrimaryVariables`

        if nargin < 5
            % Get the saturation names directly from the problem
            [~, satVars] = ...
                splitPrimaryVariables(model, problem.primaryVariables);
        end
        if isempty(satVars)
            % No saturations passed, nothing to do here.
            return
        end
        state_init = state;
        % Solution variables should be saturations directly, find the missing
        % link
        saturations = lower(model.getSaturationVarNames);
        fillsat = setdiff(saturations, lower(satVars));

        n_fill = numel(fillsat);
        assert(n_fill == 1 || n_fill == 0)
        if n_fill == 1
            fillsat = fillsat{1};
            % Fill component is whichever saturation is assumed to fill up the rest of
            % the pores. This is done by setting that increment equal to the
            % negation of all others so that sum(s) == 0 at end of update
            solvedFor = ~strcmpi(saturations, fillsat);
        else
            % All saturations are primary variables. Sum of saturations is
            % assumed to be enforced from the equation setup
            solvedFor = true(numel(saturations), 1);
        end
        ds = zeros(model.G.cells.num, numel(saturations));

        tmp = 0;
        for i = 1:numel(saturations)
            if solvedFor(i)
                v = model.getIncrement(dx, problem, saturations{i});
                ds(:, i) = v;
                if n_fill > 0
                    % Saturations added for active variables must be subtracted
                    % from the last phase
                    tmp = tmp - v;
                end
            end
        end
        ds(:, ~solvedFor) = tmp;
        % We update all saturations simultanously, since this does not bias the
        % increment towards one phase in particular.
        state = model.updateStateFromIncrement(state, ds, problem, 's', inf, model.dsMaxAbs);
        if isempty(model.FlowPropertyFunctions)
            chopped = false(model.G.cells.num, 1);
        else
            kr = model.FlowPropertyFunctions.RelativePermeability;
            [state, chopped] = kr.applyImmobileChop(model, state, state_init);
        end
        if n_fill == 1
            % Ensure that values are within zero->one interval, and
            % re-normalize if any values were capped
            bad = any((state.s > 1) | (state.s < 0) | chopped, 2);
            if any(bad)
                state.s(bad, :) = min(state.s(bad, :), 1);
                state.s(bad, :) = max(state.s(bad, :), 0);
                state.s(bad, :) = bsxfun(@rdivide, state.s(bad, :), sum(state.s(bad, :), 2));
            end
        else
            % Ensure positive values, we assume that sum of saturations is
            % handled via constraint equation
            state.s = max(state.s, 0);
        end
    end

    % --------------------------------------------------------------------%
    function state = setPhaseData(model, state, data, fld, subs)
        % Store phase data in state for further output.
        %
        % SYNOPSIS:
        %   state = model.setPhaseData(state, data, 'someField')
        %   state = model.setPhaseData(state, data, 'someField', indices)
        %
        % DESCRIPTION:
        %   Utility function for storing phase data in the state. This is
        %   used for densities, fluxes, mobilities and so on when requested
        %   from the simulator.
        %
        % PARAMETERS:
        %   model - Class instance
        %   data  - Cell array of data to be stored. One entry per active
        %           phase in `model`.
        %   fld   - The field to be stored.
        %   subs  - OPTIONAL. The subset for which phase data is to be
        %           stored. Must be a valid index of the type ::
        %             data.(fld)(subs, someIndex) = data{i}
        %           for all i. Defaults to all indices.
        %
        % RETURNS:
        %   state - state with updated `.fld`.
        %
        if nargin == 4
            subs = ':';
        end
        isActive = model.getActivePhases();

        ind = 1;
        for i = 1:numel(data)
            if isActive(i)
                state.(fld)(subs, ind) = data{i};
                ind = ind + 1;
            end
        end
    end

    % --------------------------------------------------------------------%
    function state = storeFluxes(model, state, vW, vO, vG)
        % Store integrated internal phase fluxes in state.
        %
        % SYNOPSIS:
        %   Three-phase case
        %   state = model.storeFluxes(state, vW, vO, vG);
        %   % Only water and gas in model:
        %   state = model.storeFluxes(state, vW, [], vG);
        %
        %
        % PARAMETERS:
        %   model - Class instance
        %   vW    - Water fluxes, one value per internal interface.
        %   vO    - Oil fluxes, one value per internal interface.
        %   vG    - Gas fluxes, one value per internal interface.
        %
        % RETURNS:
        %   state - State with stored values
        %
        % SEE ALSO:
        %   `storeBoundaryFluxes`

        isActive = model.getActivePhases();

        internal = model.operators.internalConn;
        state.flux = zeros(numel(internal), sum(isActive));
        phasefluxes = {value(vW), value(vO), value(vG)};
        state = model.setPhaseData(state, phasefluxes, 'flux', internal);
    end

    function state = storeBoundaryFluxes(model, state, qW, qO, qG, forces)
        % Store integrated phase fluxes on boundary in state.
        %
        % SYNOPSIS:
        %   Three-phase case
        %   state = model.storeBoundaryFluxes(state, vW, vO, vG, forces);
        %   % Only water and gas in model:
        %   state = model.storeBoundaryFluxes(state, vW, [], vG, forces);
        %
        %
        % PARAMETERS:
        %   model  - Class instance
        %   vW     - Water fluxes, one value per BC interface.
        %   vO     - Oil fluxes, one value per BC interface.
        %   vG     - Gas fluxes, one value per BC interface.
        %   forces - `drivingForces` struct containing any boundary
        %            conditions used to obtain fluxes.
        %
        % RETURNS:
        %   state - State with stored values
        %
        % SEE ALSO:
        %   `storeFluxes`

        if ~isfield(forces, 'bc') || isempty(forces.bc)
            return
        end
        phasefluxes = {value(qW), value(qO), value(qG)};
        faces = forces.bc.face;
        % Compensate for sign. Boundary fluxes have signs that correspond
        % to in/out of the reservoir. This does not necessarily correspond
        % to the neighbor structure in the grit itself.
        sgn = 1 - 2*(model.G.faces.neighbors(faces, 2) == 0);
        for i = 1:numel(phasefluxes)
            if isempty(phasefluxes{i})
                continue
            end
            phasefluxes{i} = phasefluxes{i}.*sgn;
        end
        state = model.setPhaseData(state, phasefluxes, 'flux', faces);
    end

    % --------------------------------------------------------------------%
    function state = storeMobilities(model, state, mobW, mobO, mobG)
        % Store phase mobility per-cell in state.
        %
        % SYNOPSIS:
        %   Three-phase case
        %   state = model.storebfactors(state, mobW, mobO, mobG);
        %   % Only water and gas in model:
        %   state = model.storebfactors(state, mobW, [], mobG);
        %
        %
        % PARAMETERS:
        %   model - Class instance
        %   mobW  - Water mobilities, one value per internal interface.
        %   mobO  - Oil mobilities, one value per internal interface.
        %   mobG  - Gas mobilities, one value per internal interface.
        %
        % RETURNS:
        %   state - State with stored values
        %
        % SEE ALSO:
        %   `storeMobilities`
        isActive = model.getActivePhases();

        state.mob = zeros(model.G.cells.num, sum(isActive));
        mob = {value(mobW), value(mobO), value(mobG)};
        state = model.setPhaseData(state, mob, 'mob');
    end

    % --------------------------------------------------------------------%
    function state = storeUpstreamIndices(model, state, upcw, upco, upcg)
        % Store upstream indices for each phase in state.
        %
        % SYNOPSIS:
        %   Three-phase case
        %   state = model.storeUpstreamIndices(state, upcw, upco, upcg);
        %   % Only water and gas in model:
        %   state = model.storeUpstreamIndices(state, upcw, [], upcg);
        %
        %
        % PARAMETERS:
        %   model  - Class instance
        %   upcw   - Water upwind indicator, one value per internal
        %            interface.
        %   upco   - Oil upwind indicator, one value per internal
        %            interface.
        %   upcg   - Gas upwind indicator, one value per internal
        %            interface.
        %
        % RETURNS:
        %   state - State with stored values
        %
        % SEE ALSO:
        %   `storeBoundaryFluxes`, `storeFluxes`
        isActive = model.getActivePhases();

        nInterfaces = size(model.operators.N, 1);
        state.upstreamFlag = false(nInterfaces, sum(isActive));
        mob = {upcw, upco, upcg};
        state = model.setPhaseData(state, mob, 'upstreamFlag');
    end

    % --------------------------------------------------------------------%
    function state = storeDensity(model, state, rhoW, rhoO, rhoG)
        % Store phase densities per-cell in state.
        %
        % SYNOPSIS:
        %   Three-phase case
        %   state = model.storeDensity(state, rhoW, rhoO, rhoG);
        %   % Only water and gas in model:
        %   state = model.storeDensity(state, rhoW, [], rhoG);
        %
        %
        % PARAMETERS:
        %   model - Class instance
        %   rhoW  - Water densities, one value per internal interface.
        %   rhoO  - Oil densities, one value per internal interface.
        %   rhoG  - Gas densities, one value per internal interface.
        %
        % RETURNS:
        %   state - State with stored values
        %
        % SEE ALSO:
        %   `storeMobilities`

        isActive = model.getActivePhases();

        state.rho = zeros(model.G.cells.num, sum(isActive));
        rho = {value(rhoW), value(rhoO), value(rhoG)};
        state = model.setPhaseData(state, rho, 'rho');
    end
    % --------------------------------------------------------------------%
    function state = storebfactors(model, state, bW, bO, bG)
        % Store phase reciprocal FVF per-cell in state.
        %
        % SYNOPSIS:
        %   Three-phase case
        %   state = model.storebfactors(state, bW, bO, bG);
        %   % Only water and gas in model:
        %   state = model.storebfactors(state, bW, [], bG);
        %
        %
        % PARAMETERS:
        %   model - Class instance
        %   bW  - Water reciprocal FVF, one value per internal interface.
        %   bO  - Oil reciprocal FVF, one value per internal interface.
        %   bG  - Gas reciprocal FVF, one value per internal interface.
        %
        % RETURNS:
        %   state - State with stored values
        %
        % SEE ALSO:
        %   `storeDensities`

        isActive = model.getActivePhases();

        state.bfactor = zeros(model.G.cells.num, sum(isActive));
        b = {value(bW), value(bO), value(bG)};
        state = model.setPhaseData(state, b, 'bfactor');
    end

    % --------------------------------------------------------------------%
    function index = satVarIndex(model, name)
        % Find the index of a saturation variable by name
        %
        % SYNOPSIS:
        %   index = model.satVarIndex('water')
        %
        % PARAMETERS:
        %   model - Class instance
        %   name  - Name of phase.
        %
        % RETURNS:
        %   index - Index of the phase for this model. Empty if saturation
        %           was not found.

        index = find(strcmpi(model.getSaturationVarNames, name));
    end
    
    function index = satVarIndexMatrix(model, name)
        % Find the index of a saturation variable by name
        %
        % SYNOPSIS:
        %   index = model.satVarIndex('water')
        %
        % PARAMETERS:
        %   model - Class instance
        %   name  - Name of phase.
        %
        % RETURNS:
        %   index - Index of the phase for this model. Empty if saturation
        %           was not found.

        index = find(strcmpi(model.getSaturationVarNamesMatrix, name));
    end

    % --------------------------------------------------------------------%
    function i = compVarIndex(model, name)
        % Find the index of a component variable by name
        %
        % SYNOPSIS:
        %   index = model.compVarIndex('co2')
        %
        % PARAMETERS:
        %   model - Class instance
        %   name  - Name of component.
        %
        % RETURNS:
        %   index - Index of the component for this model. Empty if
        %           saturation was not known to model.
        %
        i = find(strcmpi(model.componentVarNames, name));
    end

    % --------------------------------------------------------------------%
    function varargout = evaluateRelPerm(model, sat, varargin)
        % Evaluate relative permeability corresponding to active phases
        %
        % SYNOPSIS:
        %   % Single-phase water model
        %   krW = model.evaluateRelPerm({sW});
        %   % Two-phase oil-water model
        %   [krW, krO] = model.evaluateRelPerm({sW, sO});
        %   % Two-phase water-gas model
        %   [krW, krG] = model.evaluateRelPerm({sW, sG});
        %   % Three-phase oil-water-gas model
        %   [krW, krO, krG] = model.evaluateRelPerm({sW, sO, sG});
        %
        % PARAMETERS:
        %   model - The model
        %   sat   - Cell array containing the saturations for all active
        %           phases.
        %
        % RETURNS:
        %   varargout - One output argument per phase present,
        %               corresponding to evaluated relative permeability
        %               functions for each phase in the canonical ordering.
        %
        % SEE ALSO:
        %   `relPermWOG`, `relPermWO`, `relPermOG`, `relPermWG`
        active = model.getActivePhases();
        nph = sum(active);
        assert(nph == numel(sat), ...
        'The number of saturations must equal the number of active phases.')
        varargout = cell(1, nph);
        names = model.getPhaseNames();

        if nph > 1
            fn = ['relPerm', names];
            [varargout{:}] = model.(fn)(sat{:}, model.fluid, varargin{:});
        elseif nph == 1
            % Call fluid interface directly if single phase
            varargout{1} = model.fluid.(['kr', names])(sat{:}, varargin{:});
        end
    end

    % --------------------------------------------------------------------%
    function g = getGravityVector(model)
        % Get the gravity vector used to instantiate the model
        %
        % SYNOPSIS:
        %   g = model.getGravityVector();
        %
        %
        % PARAMETERS:
        %   model - Class instance
        %
        % RETURNS:
        %   g     - `model.G.griddim` long vector representing the gravity
        %           accleration constant along increasing depth.
        %
        % SEE ALSO:
        %   `gravity`

        if isfield(model.G, 'griddim')
            dims = 1:model.G.griddim;
        else
            dims = ':';
        end
        g = model.gravity(dims);
    end

    % --------------------------------------------------------------------%
    function gdxyz = getGravityGradient(model)
        % Get the discretized gravity contribution on faces
        %
        % SYNOPSIS:
        %   gdz = model.getGravityGradient();
        %
        %
        % PARAMETERS:
        %   model - Class instance
        %
        % RETURNS:
        %   g     - One entry of the gravity contribution per face in the
        %           grid. Does not necessarily assume that gravity is
        %           aligned with one specific direction.
        %
        % SEE ALSO:
        %   `gravity`

        assert(isfield(model.G, 'cells'), 'Missing cell field on grid');
        assert(isfield(model.G.cells, 'centroids'),...
            'Missing centroids field on grid. Consider using computeGeometry first.');

        g = model.getGravityVector();
        gdxyz = model.operators.Grad(model.G.cells.centroids) * g';
    end

% --------------------------------------------------------------------%
    function scaling = getScalingFactorsCPR(model, problem, names, solver) %#ok
        % Get scaling factors for CPR reduction in `CPRSolverAD`
        %
        % PARAMETERS:
        %   model   - Class instance
        %   problem - `LinearizedProblemAD` which is intended for CPR
        %             preconditioning.
        %   names   - The names of the equations for which the factors are
        %             to be obtained.
        %   solver  - The `LinearSolverAD` class requesting the scaling
        %             factors.
        %
        % RETURNS:
        %   scaling - Cell array with either a scalar scaling factor for
        %             each equation, or a vector of equal length to that
        %             equation.
        %
        % SEE ALSO
        %   `CPRSolverAD`
        scaling = cell(numel(names), 1);
        [scaling{:}] = deal(1);
    end

% --------------------------------------------------------------------%

    function [eqs, names, types, wellSol, src] = insertWellEquations(model, eqs, names, ...
                                                     types, wellSol0, wellSol, ...
                                                     wellVars, wellMap, ...
                                                     p, mob, rho, ...
                                                     dissolved, components, ...
                                                     dt, opt)
        % Add in the effect of wells to a system of equations, by adding
        % corresponding source terms and augmenting the system with
        % additional equations for the wells.
        %
        % PARAMETERS:
        %
        %   eqs    - Cell array of equations that are to be updated.
        %
        %   names  - The names of the equations to be updated. If
        %            phase-pseudocomponents are to be used, the names must
        %            correspond to some combination of "water", "oil", "gas"
        %            if no special component treatment is to be introduced.
        %
        %   types  - Cell array with the types of "eqs". Note that these
        %            types must be 'cell' where source terms is to be added.
        %
        %   src    - Struct containing all the different source terms that
        %            were computed and added to the equations.
        %
        %   various - Additional input arguments correspond to standard
        %             reservoir-well coupling found in `FacilityModel.
        if model.FacilityModel.getNumberOfActiveWells(wellSol) == 0
            src = [];
            return
        end
        fm = model.FacilityModel;
        nPh = nnz(model.getActivePhases);
        [src, wellsys, wellSol] = ...
            fm.getWellContributions(wellSol0, wellSol, wellVars, ...
                                    wellMap, p, mob, rho, dissolved, components, ...
                                    dt, opt.iteration);

        rhoS = model.getSurfaceDensities();
        wc = src.sourceCells;
        [~, longNames] = getPhaseNames(model);
        % Treat phase pseudocomponent source terms from wells
        for i = 1:nPh
            sub = strcmpi(names, longNames{i});
            if any(sub)
                assert(strcmpi(types{sub}, 'cell'), 'Unable to add source terms to equation that is not per cell.');
                eqs{sub}(wc) = eqs{sub}(wc) - src.phaseMass{i}./rhoS(i);
            end
        end
        % Treat component source terms from wells
        cnames = model.getComponentNames();
        for i = 1:numel(cnames)
            sub = strcmpi(names, cnames{i});
            if any(sub)
                assert(strcmpi(types{sub}, 'cell'), 'Unable to add source terms to equation that is not per cell.');
                eqs{sub}(wc) = eqs{sub}(wc) - src.components{i};
            end
        end
        eqs = horzcat(eqs, wellsys.wellEquations, {wellsys.controlEquation});
        names = horzcat(names, wellsys.names, 'closureWells');
        types = horzcat(types, wellsys.types, 'well');
    end
    function [eqs, state, src] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                                 p, s, mob, rho, ...
                                                                 dissolved, components, ...
                                                                 forces)
        % Add in the boundary conditions and source terms to equations
        %
        % SYNOPSIS:
        %   [eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
        %                                                        p, sat, mob, rho, ...
        %                                                        rs, components, ...
        %                                                        drivingForces);
        %
        % PARAMETERS:
        %   model  - Class instance.
        %   eqs    - Cell array of equations that are to be updated.
        %
        %   names  - The names of the equations to be updated. If
        %            phase-pseudocomponents are to be used, the names must
        %            correspond to some combination of "water", "oil", "gas"
        %            if no special component treatment is to be introduced.
        %
        %   types  - Cell array with the types of "eqs". Note that these
        %            types must be 'cell' where source terms is to be added.
        %
        %   src    - Struct containing all the different source terms that
        %            were computed and added to the equations.
        %
        %   p      - Cell array of phase pressures.
        %
        %   s      - Cell array of phase saturations.
        %
        %   mob    - Cell array of phase mobilities
        %
        %   rho    - Cell array of phase densities
        %
        %   dissolved - Cell array of dissolved components for black-oil
        %               style pseudocompositional models.
        %
        %   components - Cell array of equal length to the number of
        %                components. The exact representation may vary
        %                based on the model, but the respective
        %                sub-component is passed onto
        %                `addComponentContributions`.
        %
        %   forces - DrivingForces struct (see `getValidDrivingForces`)
        %            containing (possibily empty) `src` and `bc` fields.
        %
        % RETURNS:
        %   eqs   - Equations with corresponding source terms added.
        %   state - Reservoir state. Can be modified to store e.g. boundary
        %           fluxes due to boundary conditions.
        %   src   - Normalized struct containing the source terms used.
        %
        % NOTE:
        %  This function accomodates both the option of black-oil
        %  pseudocomponents (if the model equations are named "oil", "gas"
        %  or water) and true components existing in multiple phases.
        %  Mixing the two behaviors can lead to unexpected source terms.

        [src_terms, bnd_cond] = computeSourcesAndBoundaryConditionsAD(model, p, s, mob, rho, dissolved, forces);
        [~, longNames] = getPhaseNames(model);
        rhoS = model.getSurfaceDensities();
        % We first consider pseudocomponents that correspond to phases,
        % e.g. the black-oil model, most immiscible models and other models
        % where the number of phases is approximately equal to the number
        % of components. If the equations "water", "oil" and "gas" exist,
        % these will get direct source terms added. Note that the
        % corresponding source terms will already have added the effect of
        % dissolution (rs/rv).
        %
        % For fully compositional problems, this branch will not execute.
        for i = 1:numel(s)
            sub = strcmpi(names, longNames{i});
            if any(sub)
                assert(strcmpi(types{sub}, 'cell'), 'Unable to add source terms to equation that is not per cell.');
                sc = src_terms.sourceCells;
                if ~isempty(sc)
                    eqs{sub}(sc) = eqs{sub}(sc) - src_terms.phaseMass{i}./rhoS(i);
                end

                if isfield(forces.bc, 'phaseMass')
                    bnd_cond.phaseMass{i} = forces.bc.phaseMass(:, i);
                end

                bc = bnd_cond.sourceCells;
                if ~isempty(bc)
                    if isempty(bnd_cond.mapping)
                        q = bnd_cond.phaseMass{i}./rhoS(i);
                    else
                        q = (bnd_cond.mapping*bnd_cond.phaseMass{i})./rhoS(i);
                    end
                    eqs{sub}(bc) = eqs{sub}(bc) - q;
                end
            end
        end
        % Get the fluxes and store them in the state.
        if nargout > 1 && model.outputFluxes
            act = model.getActivePhases();
            tmp = cell(numel(act), 1);
            tmp(act) = bnd_cond.phaseVolume;
            state = model.storeBoundaryFluxes(state, tmp{:}, forces);
        end
        % Finally deal with actual components that exist in the different
        % phases to varying degrees.
        cnames = model.getComponentNames();
        for i = 1:numel(cnames)
            % Iterate over individual components
            name = cnames{i};
            sub = strcmpi(name, names);
            if any(sub)
                eq = eqs{sub};
            else
                eq = zeros(model.G.cells.num, 1);
            end
            C = components{i};
            assert(strcmpi(types{sub}, 'cell'), 'Unable to add source terms to equation that is not per cell.');
            % Add BC component source terms
            [eq, bnd_cond] = model.addComponentContributions(name, eq, C, bnd_cond, forces.bc);
            [eq, src_terms] = model.addComponentContributions(name, eq, C, src_terms, forces.src);
            if any(sub)
                eqs{sub} = eq;
            end
        end
        % If requested, provide the computed values for source and bc for
        % further manipulations outside this function.
        if nargout > 2
            src = struct('src', src_terms, 'bc', bnd_cond);
        end
    end

    function [eq, src] = addComponentContributions(model, cname, eq, component, src, force)
        % For a given component conservation equation, compute and add in
        % source terms for a specific source/bc where the fluxes have
        % already been computed.
        %
        % PARAMETERS:
        %
        %   model  - (Base class, automatic)
        %
        %   cname  - Name of the component. Must be a property known to the
        %            model itself through `getProp` and `getVariableField`.
        %
        %   eq     - Equation where the source terms are to be added. Should
        %            be one value per cell in the simulation grid (model.G)
        %            so that the src.sourceCells is meaningful.
        %
        %   component - Cell-wise values of the component in question. Used
        %               for outflow source terms only.
        %
        %   src    - Source struct containing fields for fluxes etc. Should
        %            be constructed from force and the current reservoir
        %            state by `computeSourcesAndBoundaryConditionsAD`.
        %
        %   force  - Force struct used to produce src. Should contain the
        %            field defining the component in question, so that the
        %            inflow of the component through the boundary condition
        %            or source terms can accurately by estimated.

            error('This function is not valid for the base class. It should be implemented for the derived classes that contain components');

    end

    function rhoS = getSurfaceDensities(model, varargin)
        % Get the surface densities of the active phases in canonical
        % ordering (WOG, with any inactive phases removed).
        %
        % RETURNS:
        %   rhoS - pvt x n double array of surface densities.
        names = model.getPhaseNames();
        rhoS = value(arrayfun(@(x) model.fluid.(['rho', x, 'S'])', names, 'UniformOutput', false));
    end

    function [compEqs, compSrc, eqNames, wellSol] = getExtraWellContributions(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration)
        % This function is called by the well model (base class:
        % SimpleWell) during the assembly of well equations and addition of
        % well source terms. The purpose of this function, given the
        % internal variables of the well model, is to compute the additional
        % closure equations and source terms that the model requires. For
        % instance, if the model contains different components that require
        % special treatment (see for example the implementation of this
        % function in :class:`ad_eor.models.OilWaterPolymerModel` in the
        % `ad-eor` module), this function should assemble any additional
        % equations and corresponding source terms. It is also possible to
        % add source terms without actually adding well equations.
        %
        % RETURNS:
        %
        %   compEqs - A cell array of additional equations added to the
        %             system to account for the treatment of components etc
        %             in the well system.
        %   compSrc - Cell array of component source terms, ordered and with
        %             the same length as the output from
        %             `ReservoirModel.getComponentNames`.
        %   eqNames - Names of the added equations. Must correspond to the
        %             same entries as `getExtraWellEquationNames` (but does
        %             not have to maintain the same ordering).
        %
        % NOTE:
        %   Input arguments are intentionally undocumented and subject to
        %   change. Please see `SimpleWell` for details.
        [compEqs, compSrc, eqNames] = deal({});
    end

    function [names, types] = getExtraWellEquationNames(model)
        % Get the names and types of extra well equations in model
        %
        % SYNOPSIS:
        %   [names, types] = model.getExtraWellEquationNames();
        %
        %
        % PARAMETERS:
        %   model - Base class.
        %
        % RETURNS:
        %   names - Cell array of additional well equations.
        %   types - Cell array of corresponding types to `names`.
        %
        % SEE ALSO:
        %   `getExtraWellContributions`

        [names, types] = deal({});
    end

    function names = getExtraWellPrimaryVariableNames(model)
        % Get the names of extra well primary variables required by model.
        %
        % SYNOPSIS:
        %   names = model.getExtraWellPrimaryVariableNames();
        %
        % PARAMETERS:
        %   model - Class instance
        %
        % RETURNS:
        %   names - Cell array of named primary variables.
        %
        % SEE ALSO:
        %   `getExtraWellContributions`
        names = {};
    end
end

methods (Static)
    % --------------------------------------------------------------------%
    function [krW, krO, krG] = relPermWOG(sw, so, sg, f, varargin)
        % Three-phase water-oil-gas relative permeability function
        %
        % SYNOPSIS:
        %   [krW, krO, krG] = model.relPermWOG(sw, so, sg, f);
        %
        %
        % PARAMETERS:
        %   sw  - Water saturation
        %   so  - Oil saturation
        %   sg  - Gas saturation
        %   f   - Struct representing the field. Fields that are used:
        %
        %           - `krW`:  Water relperm function of water saturation.
        %           - `krOW`: Oil-water relperm function of oil saturation.
        %           - `krOG`: Oil-gas relperm function of oil saturation.
        %           - `krG`:   Gas relperm function of gas saturation.
        %           - `sWcon`: Connate water saturation. OPTIONAL.
        % RETURNS:
        %   krW - Water relative permeability.
        %   krO - Oil relative permeability.
        %   krG - Gas relative permeability
        %
        % NOTE:
        %   This function should typically not be called directly as its
        %   interface is subject to change. Instead, use `evaluateRelPerm`.
        swcon = 0;
        if isfield(f, 'sWcon')
            if isempty(varargin) || numel(f.sWcon) == 1
                swcon = f.sWcon;
            else
                assert(strcmp(varargin{1}, 'cellInx'))
                swcon = f.sWcon(varargin{2});
            end
        end
        swcon = min(swcon, value(sw)-1e-5);

        d  = (sg+sw-swcon);
        ww = (sw-swcon)./d;
        krW = f.krW(sw, varargin{:});

        wg = 1-ww;
        krG = f.krG(sg, varargin{:});

        krow = f.krOW(so, varargin{:});
        krog = f.krOG(so,  varargin{:});
        krO  = wg.*krog + ww.*krow;
    end

    % --------------------------------------------------------------------%
    function [krW, krO] = relPermWO(sw, so, f, varargin)
        % Two-phase water-oil relative permeability function
        %
        % SYNOPSIS:
        %   [krW, krO] = model.relPermWO(sw, so, f);
        %
        %
        % PARAMETERS:
        %   sw  - Water saturation
        %   so  - Oil saturation
        %   f   - Struct representing the field. Fields that are used:
        %
        %            - `krW`: Water relperm function of water saturation.
        %            - `krO`: Oil relperm function of oil saturation.
        %            - `krOW`: Oil-water relperm function of oil saturation.
        %              Only used if `krO` is not found.
        % RETURNS:
        %   krW - Water relative permeability.
        %   krO - Oil relative permeability.
        %
        % NOTE:
        %   This function should typically not be called directly as its
        %   interface is subject to change. Instead, use `evaluateRelPerm`.
        krW = f.krW(sw, varargin{:});
        if isfield(f, 'krO')
            krO = f.krO(so, varargin{:});
        else
            krO = f.krOW(so, varargin{:});
        end
    end

    % --------------------------------------------------------------------%
    function [krO, krG] = relPermOG(so, sg, f, varargin)
        % Two-phase oil-gas relative permeability function
        %
        % SYNOPSIS:
        %   [krO, krG] = model.relPermOG(so, sg, f);
        %
        %
        % PARAMETERS:
        %   sw  - Water saturation
        %   sg  - Gas saturation
        %   f   - Struct representing the field. Fields that are used:
        %
        %           - `krO`: Oil relperm function of oil saturation.
        %           - `krOG`: Oil-gas relperm function of gas saturation.
        %             This function is only used if `krO` is not found.
        %           - `krG`: Gas relperm function of gas saturation.
        %
        % RETURNS:
        %   krO - Oil relative permeability.
        %   krG - Gas relative permeability
        %
        % NOTE:
        %   This function should typically not be called directly as its
        %   interface is subject to change. Instead, use `evaluateRelPerm`.
        krG = f.krG(sg, varargin{:});
        if isfield(f, 'krO')
            krO = f.krO(so, varargin{:});
        else
            krO = f.krOG(so, varargin{:});
        end
    end

    % --------------------------------------------------------------------%
    function [krW, krG] = relPermWG(sw, sg, f, varargin)
        % Two-phase water-gas relative permeability function
        %
        % SYNOPSIS:
        %   [krW, krG] = model.relPermWG(sw, sg, f);
        %
        %
        % PARAMETERS:
        %   sw  - Water saturation
        %   sg  - Gas saturation
        %   f   - Struct representing the field. Fields that are used:
        %
        %           - `krW`: Water relperm function of water saturation.
        %           - `krG`: Gas relperm function of gas saturation.
        %
        % RETURNS:
        %   krW - Water relative permeability.
        %   krG - Gas relative permeability
        %
        % NOTE:
        %   This function should typically not be called directly as its
        %   interface is subject to change. Instead, use `evaluateRelPerm`.
        krG = f.krG(sg, varargin{:});
        krW = f.krW(sw, varargin{:});
    end

    % --------------------------------------------------------------------%
    function ds = adjustStepFromSatBounds(s, ds)
        % Ensure that cellwise increment for each phase is done with
        % the same length, in a manner that avoids saturation
        % violations.
        tmp = s + ds;

        violateUpper =     max(tmp - 1, 0);
        violateLower = abs(min(tmp    , 0));

        violate = max(violateUpper, violateLower);

        [worst, jj]= max(violate, [], 2);

        bad = worst > 0;
        if any(bad)
            w = ones(size(s, 1), 1);
            for i = 1:size(s, 2)
                ind = bad & jj == i;
                dworst = abs(ds(ind, i));

                w(ind) = (dworst - worst(ind))./dworst;
            end
            ds(bad, :) = bsxfun(@times, ds(bad, :), w(bad, :));
        end
    end
end
end

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

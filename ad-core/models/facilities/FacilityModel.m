classdef FacilityModel < PhysicalModel
    % A model coupling facilities and wells to the reservoir
    %
    % SYNOPSIS:
    %
    %   model = FacilityModel(reservoirModel)
    %
    % DESCRIPTION:
    %   The `FacilityModel` is the layer between a `ReservoirModel` and the
    %   facilities. The facilities consist of a number of different wells
    %   that are implemented in their own subclasses.
    %
    %   Different wells have different governing equations, primary
    %   variables and convergence criterions. This class seamlessly handles
    %   wells appearing and disappearing, controls changing and even well
    %   type changing.
    %
    % PARAMETERS:
    %   resModel - `ReservoirModel` derived class the facilities are
    %               coupled to.
    %
    % OPTIONAL PARAMETERS:
    %   'property' - Set property to the specified value.
    %
    % RETURNS:
    %   model - Class instance of `FacilityModel`.
    %
    %
    % SEE ALSO:
    %   `ReservoirModel`, `PhysicalModel`, `SimpleWell`

    properties
        WellModels % Cell array of instantiated wells
        FacilityFlowDiscretization

        toleranceWellBHP % Convergence tolerance for BHP-type controls
        toleranceWellRate % Convergence tolerance for rate-type controls
        toleranceWellMS % Convergence tolerance for multisegment wells
        ReservoirModel % The model instance the FacilityModel is coupled to

        VFPTablesInjector % Injector VFP Tables. EXPERIMENTAL.
        VFPTablesProducer % Producer VFP Tables. EXPERIMENTAL.
        primaryVariableSet = 'standard'; % Default: bhp + phase surface rates
        setWellTargets = true; % Explicitly wellSol values to imposed targets after each iteration
    end

    properties (SetAccess = protected)
        addedPrimaryVarNames = {}; % Canonical list of all extra primary variables added by the wells
        addedPrimaryVarNamesIsFromResModel = [];  % Indicator, per primary variable, if it was added by the reservoir model (true) or if it is from the well itself (false)
        addedEquationNames = {}; % Canonical list of additional equations
        addedEquationTypes = {}; % Canonical list of the types of the added equations
    end

    methods
        function model = FacilityModel(varargin)
            model = model@PhysicalModel([]);
            if mod(nargin, 2) == 1
                model.ReservoirModel = varargin{1};
                varargin = varargin(2:end);
            end
            model.toleranceWellBHP  = 1*barsa;
            model.toleranceWellRate = 1/day;
            model.toleranceWellMS   = 1e-6;
            model.VFPTablesInjector = {};
            model.VFPTablesProducer = {};
            model = merge_options(model, varargin{:});
            model.WellModels = {};
        end
        
        function model = setupStateFunctionGroupings(model, useDefaults)
            if nargin < 2
                useDefaults = isempty(model.FacilityFlowDiscretization);
            end
            if useDefaults
                model.FacilityFlowDiscretization = FacilityFlowDiscretization(model); %#ok
            end
        end
        
        function model = removeStateFunctionGroupings(model)
            model.FacilityFlowDiscretization = [];
        end
        
        function [p, state] = getProp(facility, state, name)
            % Helper function to get (possible) control values for all
            % active wells. This function will always prioritize primary
            % variables (AD status) over doubles used for controls.
            if strcmpi(name, 'surfaceRates') || strcmpi(name, 'q_s')
                % Ask for all phase rates, bundled up
                names = facility.ReservoirModel.getPhaseNames();
                nph = numel(names);
                tmp = cell(1, nph);
                qnames = arrayfun(@(x) ['q', x, 's'], names, 'UniformOutput', false);
                [tmp{:}] = facility.getProps(state, qnames{:});
                p = tmp;
            else
                hasFState = isfield(state, 'FacilityState');
                if hasFState
                    fs = state.FacilityState;
                    pv = strcmpi(fs.names, name);
                    if any(pv)
                        % We found it!
                        p = fs.primaryVariables{pv};
                        return
                    end
                end
                % It could be stored under the variable field as well,
                % check that too
                [s, index] = facility.getVariableField(name, false);
                if hasFState
                    pv = strcmpi(fs.names, s);
                    if any(pv)
                        p = fs.primaryVariables{pv};
                        return
                    end
                end
                % We didn't find it in primary variables, fall back to
                % wellSols for non-AD values
                if isfield(state.wellSol, s)
                    wstatus = vertcat(state.wellSol.status);
                    p = arrayfun(@(x) x.(s)(:, index), state.wellSol(wstatus), 'UniformOutput', false);
                    p = vertcat(p{:});
                else
                    [p, state] = getProp@PhysicalModel(facility, state, name);
                end
            end
        end
        
        function [fn, index] = getVariableField(model, name, varargin)
            index = ':';
            names = model.ReservoirModel.getPhaseNames();
            if strcmp(name(1), 'q')
                qnames = arrayfun(@(x) ['q', x, 's'], names, 'UniformOutput', false);
            else
                qnames = {};
            end
            switch name
                case 'bhp'
                    fn = name;
                case qnames
                    fn = name;
                otherwise
                    [fn, index] = getVariableField@PhysicalModel(model, name, varargin{:});
            end
        end

        function model = setupWells(model, W, wellmodels)
            % Set up well models for changed controls or the first
            % simulation step.
            %
            % PARAMETERS:
            % 
            %   W       - Well struct (obtained from e.g. `addWell` or
            %             `processWells`)
            %
            % OPTIONAL PARAMETERS:
            %   wellmodels - Cell array of equal length to W, containing class
            %                instances for each well (e.g. `SimpleWell`,
            %                `MultisegmentWell`, or classes derived from
            %                these).  If not provided, well models be
            %                constructed from the input.
            %
            % RETURNS:
            %   model - Updated `FacilityModel` instance ready for use with
            %           wells of type `W`.
            nw = numel(W);
            if model.getNumberOfWells == 0
                % First time setup
                [pvars, fromRes, eqnames, eqtypes] = deal(cell(nw, 1));
                model.WellModels = cell(nw, 1);
                for i = 1:nw
                    if W(i).status
                        assert(numel(W(i).val) == 1, ...
                            'Bad controls for well %d: One value must be specified for .val', i);
                    end
                    % Set up models: SimpleWell or MultisegmentWell
                    if nargin < 3
                        if isfield(W(i), 'isMS') && W(i).isMS
                            wm = MultisegmentWell(W(i));
                        else
                            wm = SimpleWell(W(i));
                        end
                    else
                        wm = wellmodels{i};
                    end
                    wm.dsMaxAbs = model.ReservoirModel.dsMaxAbs;
                    wm.dpMaxAbs = model.ReservoirModel.dpMaxAbs;
                    if isfield(W(i), 'vfp_index')
                        vfp_ix = W(i).vfp_index;
                        if vfp_ix > 0
                            if wm.isInjector()
                                vfp = model.VFPTablesInjector{vfp_ix};
                            else
                                vfp = model.VFPTablesProducer{vfp_ix};
                            end
                            wm.VFPTable = vfp;
                        end
                    end
                    % Get the added primary variables for this well, plus
                    % the equations and equation types it adds
                    [pvars{i}, fromRes{i}] = wm.getExtraPrimaryVariableNames(model.ReservoirModel);
                    [eqnames{i}, eqtypes{i}] = wm.getExtraEquationNames(model.ReservoirModel);
                    model.WellModels{i} = wm;
                end
                % Combine the different equations and types added by the
                % different wells into a canonical ordering.
                [model.addedPrimaryVarNames, keepix] = uniqueStable([pvars{:}]);
                model.addedPrimaryVarNamesIsFromResModel = [fromRes{:}];
                model.addedPrimaryVarNamesIsFromResModel = model.addedPrimaryVarNamesIsFromResModel(keepix);
                
                [model.addedEquationNames, keepix] = uniqueStable([eqnames{:}]);

                etypes = [eqtypes{:}];
                model.addedEquationTypes = etypes(keepix);
            else
                assert(model.getNumberOfWells == nw, ...
                    'Number of wells in facility model has changed during simulation')
                for i = 1:nw
                    % Update with new wells. Typically just a quick
                    % overwrite of existing wells
                    model.WellModels{i} = model.WellModels{i}.updateWell(W(i));
                end
            end
        end

        function model = setReservoirModel(model, resmodel)
            model.ReservoirModel = resmodel;
            model.AutoDiffBackend = resmodel.AutoDiffBackend;
        end
        function W = getWellStruct(model, subs)
            % Get the well struct representing the current set of wells
            %
            % SYNOPSIS:
            %   W = model.getWellStruct();
            %
            % PARAMETERS:
            %   model - `FacilityModel` class instance
            %
            % RETURNS:
            %   W - Standard well struct.
            %
            if nargin == 1
                subs = ':';
            end
            W = cellfun(@(x) x.W, model.WellModels(subs), 'UniformOutput', false);
            W = vertcat(W{:});
        end
        
        function nwell = getNumberOfActiveWells(model, wellSol)
            % Get number of wells active initialized in facility
            %
            % SYNOPSIS:
            %   W = model.getNumberOfActiveWells();
            %
            % PARAMETERS:
            %   model - `FacilityModel` class instance
            %
            % RETURNS:
            %   nwell - Number of active wells
            %
            mask = model.getWellStatusMask(wellSol);
            nwell = nnz(mask);
        end

        function nwell = getNumberOfWells(model)
            % Get number of wells initialized in facility
            %
            % SYNOPSIS:
            %   W = model.getNumberOfWells();
            %
            % PARAMETERS:
            %   model - `FacilityModel` class instance
            %
            % RETURNS:
            %   nwell - Number of wells
            %
            nwell = numel(model.WellModels);
        end
        
        function actIx = getIndicesOfActiveWells(model, wellSol)
            % Get indices of active wells
            %
            % SYNOPSIS:
            %   actIx = model.getIndicesOfActiveWells(wellSol);
            %
            % PARAMETERS:
            %   model   - `FacilityModel` class instance
            %   wellSol - The wellSol struct
            %
            % RETURNS:
            %   actIx - The indices of the active wells in the global list
            %           of all wells (active & inactive)
            %
            act = model.getWellStatusMask(wellSol);
            actIx = (1:numel(act))';
            actIx = actIx(act);
        end

        function act = getWellStatusMask(model, wellSol)
            % Get status mask for active wells
            %
            % SYNOPSIS:
            %   act = model.getWellStatusMask(wellSol);
            %
            % DESCRIPTION:
            %   Get the well status of all wells. The status is true if the
            %   well is present and active. Wells can be disabled in two
            %   ways: Their status flag can be set to false in the well
            %   struct, or the wellSol.status flag can be set to false by 
            %   the simulator itself.
            %
            % PARAMETERS:
            %   model   - `FacilityModel` class instance
            %   wellSol - The wellSol struct
            %
            % RETURNS:
            %   act   - Array with equal length to the total number of
            %           wells, with booleans indicating if a specific well
            %           is currently active.
            %
            actModel = cellfun(@(x) x.W.status, model.WellModels);
            actWellSol = arrayfun(@(x) x.status, wellSol);
            
            act = reshape(actModel, [], 1) & reshape(actWellSol, [], 1);
        end

        function state = initStateAD(model, state, vars, names, origin)
            isComp = strcmp(names, 'well_massfractions');
            if any(isComp)
                mf = vars{isComp};
                vars = vars(~isComp);
                names = names(~isComp);
            end
            phnames = model.ReservoirModel.getPhaseNames();
            nph = numel(phnames);
            surfaceRates = cell(1, nph);
            act = vertcat(state.wellSol.status);
            for i = 1:nph
                fn = ['q', phnames(i), 's'];
                sub = strcmp(names, fn);
                if any(sub)
                    surfaceRates{i} = vars{sub};
                else
                    surfaceRates{i} = vertcat(state.wellSol(act).(fn));
                end
            end
            state.FacilityState = struct('primaryVariables', {vars}, ...
                                         'surfacePhaseRates', {surfaceRates}, ...
                                         'names', {names});
            if any(isComp)
                ncomp = model.ReservoirModel.getNumberOfComponents();
                state.FacilityState.massfractions = cell(1, ncomp);
                len = numelValue(mf)/(ncomp-1);
                fill = 1;
                for i = 0:ncomp-2
                    index = len*i + (1:len);
                    v = mf(index);
                    fill = fill - v;
                    state.FacilityState.massfractions{i+1} = v;
                end
                state.FacilityState.massfractions{end} = fill;
            end
            state = initStateAD@PhysicalModel(model, state, {}, {}, {});
        end

        function [variables, names, origin] = getPrimaryVariables(model, state)
            if isfield(state, 'wellSol')
                [variables, names] = model.getBasicPrimaryVariables(state.wellSol);
                origin = cell(size(variables));
                [origin{:}] = deal(class(model));
            else
                [variables, names, origin] = deal({});
            end
        end

        function model = validateModel(model, varargin)
            if nargin > 1
                W = varargin{1}.W;
                model = model.setupWells(W);
            end
        end
        
        function names = getPrimaryVariableNames(model)
            % Get the names of primary variables present in all wells
            %
            % SYNOPSIS:
            %   names = model.getPrimaryVariableNames();
            %
            % DESCRIPTION:
            %   Get a list of the names of primary variables that will be
            %   required to solve a problem with the current set of wells.
            %
            % PARAMETERS:
            %   model  - `FacilityModel` class instance
            %
            % RETURNS:
            %   names  - Cell array of the names of the primary variables.
            %
            names = [model.getBasicPrimaryVariableNames(), model.addedPrimaryVarNames];
        end

        function names = getBasicPrimaryVariableNames(model)
            % Get the names of the basic primary variables present in all wells
            %
            % SYNOPSIS:
            %   names = model.getBasicPrimaryVariableNames();
            %
            % DESCRIPTION:
            %   Get a list of the basic names of primary variables that will 
            %   be required to solve a problem with the current set of wells.
            %   The basic primary variables are always present in MRST, and
            %   correspond to the phase rates for each phase present, as
            %   well as the bottom-hole pressures. This ensures that all
            %   solvers have a minimum feature set for well controls.
            %
            % PARAMETERS:
            %   model  - `FacilityModel` class instance
            %
            % RETURNS:
            %   names  - Cell array of the names of the basic primary
            %            variables.
            %
            switch lower(model.primaryVariableSet)
                case {'standard', 'explicit'}
                    phNames = model.ReservoirModel.getPhaseNames();
                    names = arrayfun(@(x) ['q', x, 's'], phNames, 'UniformOutput', false);
                    names = [names, 'bhp'];
                case 'bhp'
                    names = {'bhp'};
                case 'bhp_massfractions'
                    names = {'bhp', 'well_massfractions'};
                case 'none'
                    names = {};
                otherwise
                    error('Invalid primary variable set %s for FacilityModel',...
                                                        model.primaryVariableSet);
            end
        end

        function [variables, names, map] = getBasicPrimaryVariables(model, wellSol)
            % Get the basic primary variables common to all well models.
            %
            % SYNOPSIS:
            %   [vars, names, map] = model.getBasicPrimaryVariables(wellSol)
            %
            % DESCRIPTION:
            %   Get phase rates for the active phases and the bhp.
            %   In addition, the map contains indicators used to
            %   find the phase rates and BHP values in "variables"
            %   since these are of special importance to many
            %   applications and are considered canonical (i.e. they
            %   are always solution variables in MRST, and functions
            %   can assume that they will always be found in the
            %   variable set for wells).
            %
            % PARAMETERS:
            %   model - `FacilityModel` class instance
            %
            % RETURNS:
            %   variables - Cell array of the primary variables.
            %   names     - Cell array with the names of the primary
            %               variables.
            %   map       - Struct with details on which variables
            %               correspond to ordered phase rates and the
            %               bottom hole pressures.
            %

            
            if model.getNumberOfActiveWells(wellSol) == 0
                [variables, names] = deal({});
                [isBHP, isRate] = deal([]);
            else
                % Take the active wellSols
                active = model.getWellStatusMask(wellSol);
                wellSol = wellSol(active);
                names = model.getBasicPrimaryVariableNames();
                nvar = numel(names);
                variables = cell(1, nvar);
                isCOMP = false(size(variables));
                for i = 1:numel(variables)
                    if strcmp(names{i}, 'well_massfractions')
                        comp = vertcat(wellSol.massfractions);
                        isCOMP(i) = true;
                        variables{i} = reshape(comp(:, 1:end-1), [], 1);
                    elseif strcmp(names{i}, 'well_temperature')
                        variables{i} = vertcat(wellSol.T);
                    else
                        variables{i} = vertcat(wellSol.(names{i}));
                    end
                end
                
                isBHP = false(1, nvar);
                if nvar > 0
                    isBHP(end) = true;
                end
                isRate = ~isBHP;
            end
            map = struct('isBHP', isBHP, 'isRate', isRate);
        end

        function [variables, names, wellmap] = getExtraPrimaryVariables(model, wellSol)
            % Get additional primary variables (not in the basic set)
            %
            % SYNOPSIS:
            %   [variables, names, wellmap] = model.getExtraPrimaryVariables(wellSol)
            %
            % DESCRIPTION:
            %   Get extra primary variables are variables required by more
            %   advanced wells that are in addition to the basic facility
            %   variables (rates + bhp).
            %
            % PARAMETERS:
            %   model   - `FacilityModel` class instance
            %   wellSol - The wellSol struct
            %
            %
            % RETURNS:
            %   names     - Column of cells. Each cell is a string with the name
            %               of an extra-variable.
            %
            %   variables - Column of cells. Each element, variables{i}, is a vector given
            %               the value corresponding to extra-variable with name
            %               names{i}. This vector is composed of stacked values
            %               over all the wells that contains this extra-variable.
            %
            %   wellmap   - The facility model contains the extra-variables of all
            %               the well models that are used. Let us consider the well
            %               with well number wno (in the set of active wells), then
            %               the Well model is belongs to has its own
            %               extra-variables (a subset of those of the Facility
            %               model). We consider the j-th extra-variable of
            %               the Well model. Then, `i = extraMap(wno, j)`
            %               says that this extra-variable corresponds to
            %               `names{i}`.
            %
            % SEE ALSO:
            %   `getBasicPrimaryVariables`
            actIx = model.getIndicesOfActiveWells(wellSol);
            nw = numel(actIx);
            if nw == 0
                [variables, names, wellmap] = deal({});
                return
            end
            names = model.addedPrimaryVarNames;
            nv = numel(names);
            vars = cell(nw, nv);

            wellmap = zeros(nw, nv);
            if nv > 0
                all_ix = (1:nv)';
                for i = 1:nw
                    wno = actIx(i);
                    [v, n] = model.WellModels{wno}.getExtraPrimaryVariables(wellSol(wno), model.ReservoirModel);

                    for j = 1 : numel(v)
                        % Map into array of added primary variables
                        ix = strcmpi(names, n{j});
                        wellmap(i, j) = all_ix(ix);
                        vars{i, ix} = v{j};
                    end
                end
            end

            variables = cell(1, nv);
            for j = 1:nv
                variables{j} = vertcat(vars{:, j});
            end
        end

        function [variables, names, wellMap] = getAllPrimaryVariables(model, wellSol)
            % Get all primary variables (basic + extra)
            %
            % SYNOPSIS:
            %   [variables, names, wellmap] = model.getAllPrimaryVariables(wellSol)
            %
            % DESCRIPTION:
            %   Gets all primary variables, both basic (rates, bhp) and added
            %   variables (added by different wells and from the model
            %   itself).
            %
            % PARAMETERS:
            %   model   - `FacilityModel` class instance
            %   wellSol - The wellSol struct
            %
            % RETURNS:
            %   names     - Cell array. Each cell is a string with the name
            %               of an variable.
            %
            %   variables - Column of cells. Each element, variables{i}, is a vector given
            %               the value corresponding to variable with name
            %               names{i}. This vector is composed of stacked values
            %               over all the wells that contains this variable.
            %
            %   wellmap   - A combined struct containing mappings for both
            %               the standard and extra primary variables.
            %
            % SEE ALSO:
            %   `getBasicPrimaryVariables`, `getExtraPrimaryVariables`



            [basic, bnames, wellMap] = model.getBasicPrimaryVariables(wellSol);
            [extra, enames, extraMap] = model.getExtraPrimaryVariables(wellSol);
            
            wellMap.isBHP = [wellMap.isBHP, false(size(extra))];
            wellMap.isRate = [wellMap.isRate, false(size(extra))];
            
            wellMap.extraMap = extraMap;
            names = [bnames, enames];
            variables = [basic, extra];
        end

        function [sources, wellSystem, wellSol] = getWellContributions(model, wellSol0, wellSol, wellvars, wellMap, p, mob, rho, dissolved, comp, dt, iteration)
            % Get sources, well equations and updated wellSol
            %
            % SYNOPSIS:
            %   [sources, wellSystem, wellSol] = fm.getWellContributions(...
            %   wellSol0, wellSol, wellvars, wellMap, p, mob, rho, dissolved, comp, dt, iteration)
            %
            % DESCRIPTION:
            %   Get the source terms due to the wells, control and well
            %   equations and updated well sol. Main gateway for adding wells
            %   to a set of equations.
            %
            % PARAMETERS:
            %   model     - Facility model class instance.
            %   wellSol0  - wellSol struct at previous time-step.
            %   wellSol   - wellSol struct at current time-step.
            %   wellvars  - Well variables. Output from
            %               `getAllPrimaryVariables`.
            %   wellMap   - Well mapping. Output from
            %               `getAllPrimaryVariables`.
            %   p         - Pressure defined in all cells of the underlying
            %               `ReservoirModel`. Normally, this is the oil
            %               pressure.
            %   mob       - Cell array of phase mobilities defined in all
            %               cells of the reservoir.
            %   rho       - Cell array of phase densities defined in all
            %               cells of the reservoir.
            %   dissolved - Black-oil style dissolution. See
            %               :meth:`ad_blackoil.models.ThreePhaseBlackoilModel.getDissolutionMatrix`.
            %   comp      - Cell array of components in the reservoir.
            %   dt        - The time-step.
            %   iteration - The current nonlinear iteration for which the
            %               sources are to be computed.
            %
            % RETURNS:
            %   sources    - Struct containing source terms for phases,
            %                components and the corresponding cells.
            %   wellSystem - Struct containing the well equations
            %                (reservoir to wellbore, and
            %                control-equations).
            %   wellSol    - Updated wellSol struct.
            %

            if isnan(iteration) || iteration < 0
                warning(['Iteration number is not passed on to well model,', ...
                         'this may indicate wellbore pressure-drop will never be updated']);
            end
            actWellIx = model.getIndicesOfActiveWells(wellSol);
            nw = numel(actWellIx);

            % Stores base well equations describing reservoir coupling
            allBaseEqs = cell(nw, 1);
            % Control equations ensure that we enforce constraints
            allCtrl = cell(nw, 1);
            % Volumetric phase source terms
            allVol = cell(nw, 1);
            % Mass phase source terms
            allMass = cell(nw, 1);
            % Composition source terms
            allComp = cell(nw, 1);
            

            % Get the additional equations not implemented in the minimal
            % "SimpleWell" class.
            enames = model.addedEquationNames;
            etypes = model.addedEquationTypes;
            cnames = model.ReservoirModel.getComponentNames();
            ncomp = numel(cnames);
            assert(ncomp == numel(comp), ...
                'Number of input components must match length of getComponentNames!');

            n_extra = numel(enames);
            assert(numel(etypes) == n_extra);

            allExtraEqs = cell(nw, n_extra);
            resModel = model.ReservoirModel;

            addedVars = model.addedPrimaryVarNames;
            varmaps = cell(1, numel(addedVars));
            for varNo = 1:numel(addedVars)
                varmaps{varNo} = model.getWellVariableMap(addedVars{varNo}, wellSol);
            end
            
            isBH = wellMap.isBHP;
            isQ = wellMap.isRate;
            emap = wellMap.extraMap;
            
            bhp = wellvars{isBH};
            qWell = wellvars(isQ);
            wellvars = wellvars(~(isBH | isQ));

            [basenames, basetypes] = model.WellModels{1}.getWellEquationNames(resModel);
            for i = 1:nw
                wellNo = actWellIx(i);
                wm = model.WellModels{wellNo};
                ws = wellSol(wellNo);
                ws0 = wellSol0(wellNo);

                W = wm.W;
                packed = packPerforationProperties(W, p, mob, rho, dissolved, comp, wellvars, addedVars, varmaps, emap, i);
                qw = cellfun(@(x) x(i), qWell, 'uniformoutput', false);
                bh = bhp(i);
                % Update pressure
                ws = wm.updateConnectionPressureDrop(ws0, ws, resModel, qw, bh, packed, dt, iteration);
                % Update limits
                [qw, bh, ws, ok] = wm.updateLimits(ws0, ws, resModel, qw, bh, packed, dt, iteration);
                if ~ok
                    bhp(i) = bh;
                    for phNo = 1:numel(qw)
                        qWell{phNo}(i) = qw{phNo};
                    end
                end
               % Set up well equations and phase source terms
               [allBaseEqs{i}, allCtrl{i}, extraEqs, extraNames, allMass{i}, allVol{i}, ws] =...
                   wm.computeWellEquations(ws0, ws, resModel, qw, bh, packed, dt, iteration);

               % Get component source terms and corresponding equations (if
               % any components are present)
               [compEqs, allComp{i}, compNames, ws] =...
                   wm.computeComponentContributions(ws0, ws, resModel, qw, bh, packed, allMass{i}, allVol{i}, dt, iteration);

               extraEqs = {extraEqs{:}, compEqs{:}};
               extraNames = {extraNames{:}, compNames{:}};

               for eqNo = 1:numel(extraEqs)
                   % Map into global list of equations
                   ix = strcmpi(enames, extraNames{eqNo});
                   allExtraEqs{i, ix} = extraEqs{eqNo};
               end
               wellSol(wellNo) = ws;
            end
            % We have assembled all equations for each well. Combine the
            % equations from the different wells into one (array) of each
            % type.
            nPh = nnz(resModel.getActivePhases);
            [srcMass, srcVol, eqs] = deal(cell(1, nPh));
            for phNo = 1:nPh
                srcMass{phNo} = model.combineCellData(allMass, phNo);
                srcVol{phNo} = model.combineCellData(allVol, phNo);
                eqs{phNo} = model.combineCellData(allBaseEqs, phNo);
            end
            % Components are ordered canonically by reservoir model
            srcComp = cell(1, ncomp);
            for cNo = 1:ncomp
                srcComp{cNo} = model.combineCellData(allComp, cNo);
            end
            % If we have extra equations, add them in
            if n_extra > 0
                extraEqs = cell(1, n_extra);
                for i = 1:n_extra
                    ok = ~cellfun(@isempty, allExtraEqs(:, i));
                    extraEqs{i} = vertcat(allExtraEqs{ok, i});
                end
                % Equations are the base, common variables as well as any extra
                % equations added due to complex wells.
                names = horzcat(basenames, enames);
                types = horzcat(basetypes, etypes);

                eqs = {eqs{:}, extraEqs{:}};
            else
                names = basenames;
                types = basetypes;
            end
            ctrleq = vertcat(allCtrl{:});

            wc = model.getWellCells(actWellIx);
            [wc, srcMass, srcVol] = model.handleRepeatedPerforatedcells(wc, srcMass, srcVol);
            wellSystem = struct('wellEquations', {eqs}, ...
                                'names',  {names}, ...
                                'types', {types}, ...
                                'controlEquation', ctrleq);
            sources = struct('phaseMass',   {srcMass}, ...
                             'phaseVolume', {srcVol}, ...
                             'components',  {srcComp}, ...
                             'sourceCells', wc);
            if model.ReservoirModel.extraWellSolOutput
                wellSol = model.setWellSolStatistics(wellSol, sources);
            end
        end

        function wellSol = updateWellSolAfterStep(model, wellSol, wellSol0)
            % Update wellSol after step (check for closed wells, etc)
            %
            % SYNOPSIS:
            %   wellSol = model.updateWellSolAfterStep(wellSol, wellSol0)
            %
            % PARAMETERS:
            %   model     - Facility model class instance.
            %   wellSol0  - wellSol struct at previous time-step.
            %   wellSol   - wellSol struct at current time-step.
            %
            % RETURNS:
            %   wellSol   - Updated wellSol struct.
            %
            for wno = 1:numel(wellSol)
                wm = model.WellModels{wno};
                wellSol(wno) = wm.updateWellSolAfterStep(model.ReservoirModel, wellSol(wno), wellSol0(wno));
            end
        end
        
    function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
        % Generic update function for reservoir models containing wells.
        %
        % SEE ALSO:
        %   :meth:`ad_core.models.PhysicalModel.updateAfterConvergence`

        state.wellSol = model.updateWellSolAfterStep(state.wellSol, state0.wellSol);
        report = [];
    end

        function wc = getWellCells(model, subs)
            % Get the perforated cells of all wells, regardless of status
            %
            % SYNOPSIS:
            %   wc = model.getWellCells()
            %
            % PARAMETERS:
            %   model     - Facility model class instance.
            %
            % RETURNS:
            %   wc   - Array of well cells
            %
            if nargin == 1
                c = cellfun(@(x) x.W.cells, model.WellModels, 'UniformOutput', false);
            else
                c = cellfun(@(x) x.W.cells, model.WellModels(subs), 'UniformOutput', false);
            end
            wc = vertcat(c{:});
        end

        function [wc, p2w] = getActiveWellCells(model, wellSol)
            % Get the perforated cells in active wells
            %
            % SYNOPSIS:
            %   wc = model.getActiveWellCells()
            %
            % PARAMETERS:
            %   model - Facility model class instance.
            %
            % RETURNS:
            %   wc   - Array of well cells that are active, and belong to
            %          active wells.
            %
            c = cellfun(@(x) x.W.cells, model.WellModels, 'UniformOutput', false);
            active = model.getWellStatusMask(wellSol);
            wc = vertcat(c{active});
            if nargout > 1
                p2w = getPerforationToWellMapping(model.getWellStruct());
                % Map into active wells
                p2w = p2w(active(p2w));
            end
        end

        function ws = setWellSolStatistics(model, ws, sources)
            % Add statistics to wellSol (wcut, gor, ...)
            %
            % SYNOPSIS:
            %   wellSol = model.setWellSolStatistics(wellSol, sources)
            %
            % PARAMETERS:
            %   model     - Facility model class instance.
            %   wellSol   - wellSol struct at current time-step.
            %   sources   - Source struct from `getWellContributions`.
            %
            % RETURNS:
            %   wellSol   - Updated wellSol struct where additional useful
            %               information has been added
            %
            p2w = getPerforationToWellMapping(model.getWellStruct());
            % Map into active wells
            active = model.getWellStatusMask(ws);
            p2w = p2w(active(p2w));

            gind = model.ReservoirModel.getPhaseIndex('G');
            oind = model.ReservoirModel.getPhaseIndex('O');
            wind = model.ReservoirModel.getPhaseIndex('W');
            srcRes = sources.phaseVolume;
            for i = 1:numel(srcRes)
                srcRes{i} = value(srcRes{i});
            end
            qR = [srcRes{:}];
            if size(qR, 1) ~= numel(p2w)
                % Multiple wells perforated in same block, etc. Output from
                % this function is not meaningful.
                return
            end
            for i = 1:numel(ws)
                if ~active(i)
                    continue
                end
                % Store reservoir fluxes and total fluxes
                ws(i).qTs = 0;
                ws(i).qTr = 0;
                if model.ReservoirModel.gas
                    tmp = sum(srcRes{gind}(p2w == i));
                    ws(i).qGr = tmp;
                    ws(i).qTr = ws(i).qTr + tmp;
                    ws(i).qTs = ws(i).qTs + ws(i).qGs;
                end

                if model.ReservoirModel.oil
                    tmp = sum(srcRes{oind}(p2w == i));
                    ws(i).qOr = tmp;
                    ws(i).qTr = ws(i).qTr + tmp;
                    ws(i).qTs = ws(i).qTs + ws(i).qOs;
                end

                if model.ReservoirModel.water
                    tmp = sum(srcRes{wind}(p2w == i));
                    ws(i).qWr = tmp;
                    ws(i).qTr = ws(i).qTr + tmp;
                    ws(i).qTs = ws(i).qTs + ws(i).qWs;
                end

                % Phase cuts - fraction of reservoir conditions
                if model.ReservoirModel.water && model.ReservoirModel.oil
                    ws(i).wcut = ws(i).qWs./(ws(i).qWs + ws(i).qOs);
                end

                if model.ReservoirModel.gas
                    ws(i).gcut = ws(i).qGs./ws(i).qTs;
                end

                if model.ReservoirModel.oil
                    ws(i).ocut = ws(i).qOs./ws(i).qTs;
                end

                % Gas/oil ratio
                if model.ReservoirModel.gas && model.ReservoirModel.oil
                    ws(i).gor = ws(i).qGs/ws(i).qOs;
                end

                ws(i).flux = qR(p2w == i, :);
            end
        end

        % Implementation details for stand-alone model
        function [wellSol, restVars] = updateWellSol(model, wellSol, problem, dx, drivingForces, restVars)
            % Update the wellSol based on increments
            %
            % SYNOPSIS:
            %   [wellSol, restVars] = model.updateWellSol(wellSol, problem, dx, forces, restVars)
            %
            % PARAMETERS:
            %   model    - Facility model class instance
            %   wellSol  - Well solution struct
            %   problem  - Linearized problem used to produce dx.
            %   dx       - Increments corresponding to
            %             `problem.primaryVariables`
            %   forces   - Boundary condition struct
            %   restVars - Variables that have not yet been updated.
            %
            % RETURNS:
            %   state    - Updated Well solution struct
            %   restVars - Variables that have not yet been updated.

            if nargin < 6
                restVars = problem.primaryVariables;
            end
            % Update the wellSol struct
            if numel(wellSol) == 0
                % Nothing to be done if there are no wells
                return
            end
            wellVars = model.getPrimaryVariableNames();

            nVars = numel(wellVars);
            actIx = model.getIndicesOfActiveWells(wellSol);
            nW = numel(actIx);
            activeVars = false(nW, nVars);
            dx_well = cell(nW, nVars);
            % Loop over all possible well variables
            for varNo = 1:nVars
                wf = wellVars{varNo};
                if strcmp(wf, 'well_massfractions')
                    % Special case
                    c = dx{varNo};
                    c = reshape(c, nW, []);
                    c = [c, -sum(c, 2)];
                    activeVars(:, varNo) = true;
                    restVars = model.stripVars(restVars, wf);
                    for i = 1:nW
                        dx_well{i, varNo} = c(i, :);
                    end
                elseif any(strcmp(restVars, wf))
                    dv = model.getIncrement(dx, problem, wf);

                    isVarWell = model.getWellVariableMap(wf, wellSol);
                    for i = 1:nW
                        % Check if this well has this specific variable as
                        % a primary variable.
                        subs = isVarWell == i;
                        if any(subs)
                            activeVars(i, varNo) = true;
                            dx_well{i, varNo} = dv(subs);
                        end
                    end
                    % Field is taken care of
                    restVars = model.stripVars(restVars, wf);
                end
            end
            for i = 1:nW
                % Finally, update the actual wellSols using the extracted
                % increments per well.
                wNo = actIx(i);
                act = activeVars(i, :);
                dxw = dx_well(i, act);
                wellSol(wNo) = model.WellModels{wNo}.updateWellSol(wellSol(wNo), wellVars(act), dxw, model.ReservoirModel);
            end
        end

        function isVarWell = getWellVariableMap(model, wf, ws)
            % Get mapping indicating which variable belong to each well
            %
            % SYNOPSIS:
            %   isVarWell = model.getWellVariableMap('someVar', wellSol)
            %
            %
            % PARAMETERS:
            %   model - Class instance of `FacilityModel`
            %   wf    - String of variable for which the mapping will be
            %           generated.
            %   ws    - Current wellSol.
            %
            % RETURNS:
            %   isVarWell - Array equal in length to the total number of
            %               variables with name `wf`. The entries
            %               correspond to which well owns that specific
            %               variable number. This allows multiple wells to
            %               have for example bottom-hole pressures as a
            %               variable, without having to split them up by
            %               name in the reservoir equations.

            act = model.getWellStatusMask(ws);
            isRes = model.addedPrimaryVarNamesIsFromResModel;
            if isRes(strcmpi(wf, model.addedPrimaryVarNames))
                % To be expanded.
                counts = ones(size(act));
            else
                counts = cellfun(@(x) x.getVariableCounts(wf), model.WellModels(act));
            end
            isVarWell = rldecode((1:nnz(act))', counts);
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            % Get stand-alone equations for the wells
            %
            % SYNOPSIS:
            %   [problem, state] = model.getEquations(state0, state, dt, drivingForces, varargin)
            %
            % DESCRIPTION:
            %   The well equations can be solved as a separate nonlinear
            %   system with the reservoir as a fixed quantity. This is
            %   useful for debugging.
            %
            % PARAMETERS:
            %   model     - `FacilityModel` instance.
            %   wellSol0  - wellSol struct at previous time-step.
            %   wellSol   - wellSol struct at current time-step.
            %   dt        - Time-step.
            %   forces    - Forces struct for the wells.
            %
            % OPTIONAL PARAMETERS:
            %   'resOnly'  -  If supported by the equations, this flag will
            %                 result in only the values of the equations being
            %                 computed, omitting any Jacobians.
            %
            %   'iteration' - The nonlinear iteration number. This can be
            %                 provided so that the underlying equations can
            %                 account for the progress of the nonlinear
            %                 solution process in a limited degree, for example
            %                 by updating some quantities only at the first
            %                 iteration.
            % RETURNS:
            %   problem - Instance of the wrapper class `LinearizedProblemAD`
            %             containing the residual equations as well as
            %             other useful information.
            %
            %   state   - The equations are allowed to modify the system
            %             state, allowing a limited caching of expensive
            %             calculations only performed when necessary.


            opt = struct('iteration', nan, 'resOnly', false);
            opt = merge_options(opt, varargin{:});
            wellSol = state.wellSol;
            wellSol0 = state0.wellSol;
            resmodel = model.ReservoirModel;
            % Get variables from facility and wells
            [wellVars, primaryVars, wellMap] = model.getAllPrimaryVariables(wellSol);
            if ~opt.resOnly
                [wellVars{:}] = initVariablesADI(wellVars{:});
            end

            if isa(resmodel, 'ThreePhaseBlackOilModel')
                [rs, rv] = resmodel.getProps(state, 'rs', 'rv');
            else
                [rs, rv] = deal([]);
            end

            if ~isfield(state, 'rho') || ~isfield(state, 'mob')
                resmodel.FacilityModel = model;
                resmodel.extraStateOutput = true;
                % Trust that the base reservoir model has implemented
                % storage of extra properties (mobility and rho included).
                [~, state] = resmodel.getEquations(state0, state, dt, drivingForces, 'iteration', inf, 'resOnly', true);
                assert(isfield(state, 'rho'), 'Density missing from state!');
                assert(isfield(state, 'mob'), 'Mobility missing from state!');
            end
            p = resmodel.getProp(state, 'pressure');

            nPh = size(state.rho, 2);
            [mob, rho] = deal(cell(1, nPh));
            for i = 1:nPh
                mob{i} = state.mob(:, i);
                rho{i} = state.rho(:, i);
            end
            dissolution = resmodel.getDissolutionMatrix(rs, rv);
            % Note! Currently not valid for polymer or compositional
            components = {};
            [src, wellsys, state.wellSol] = ...
                model.getWellContributions(wellSol0, wellSol, wellVars, ...
                        wellMap, p, mob, rho, dissolution, components, dt, opt.iteration);

            eqs = {wellsys.wellEquations{:}, wellsys.controlEquation};
            names = {wellsys.names{:}, 'closureWells'};
            types = {wellsys.types{:}, 'well'};

            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
        end

        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            % Update state.
            %
            % SEE ALSO:
            %   :meth:`ad_core.models.PhysicalModel.updateState`

            if isfield(state, 'wellSol')
                state.wellSol = model.updateWellSol(state.wellSol, problem, dx, drivingForces);
                % Handle the directly assigned values (i.e. can be deduced directly from
                % the well controls.
                if model.setWellTargets
                    W = drivingForces.W;
                    phIndices = model.ReservoirModel.getPhaseIndices();
                    state.wellSol = assignWellValuesFromControl(model.ReservoirModel, state.wellSol, W, phIndices(1), phIndices(2), phIndices(3));
                end
            end
            report = [];
        end

        function state = validateState(model, state)
            % Validate state.
            %
            % SEE ALSO:
            %   :meth:`ad_core.models.PhysicalModel.validateState`
            if ~isfield(state, 'wellSol') || isempty(state.wellSol) || ~isfield(state.wellSol, 'bhp')
                if isfield(state, 'wellSol')
                    state = rmfield(state, 'wellSol');
                end
                W = model.getWellStruct();
                state.wellSol = initWellSolAD(W, model.ReservoirModel, state);
            end
            
            if strcmpi(model.primaryVariableSet, 'bhp_massfractions')
                ncomp = model.ReservoirModel.getNumberOfComponents();
                for i = 1:numel(state.wellSol)
                    state.wellSol(i).massfractions = ones(1, ncomp)/ncomp;
                end
            end

            for wno = 1:numel(model.WellModels)
                new_ws = model.WellModels{wno}.validateWellSol(model.ReservoirModel, state.wellSol(wno), state);
                % Hack to avoid adding fields manually
                flds = fieldnames(new_ws);
                for j = 1:numel(flds)
                    state.wellSol(wno).(flds{j}) = new_ws.(flds{j});
                end
            end
        end

        function forces = getValidDrivingForces(model)
            % Get valid driving forces.
            %
            % SEE ALSO:
            %   :meth:`ad_core.models.PhysicalModel.getValidDrivingForces`

            forces = getValidDrivingForces@PhysicalModel(model);
            forces.W   = [];
            forces.bc  = [];
            forces.src = [];
        end

        function [convergence, values, names, evaluated] = checkFacilityConvergence(model, problem)
            % Check if facility and submodels has converged
            %
            % SYNOPSIS:
            %   [convergence, values, names, evaluated] = model. checkFacilityConvergence(problem)
            %
            %
            % PARAMETERS:
            %   model   - Class instance.
            %   problem - `LinearizedProblemAD` to be checked for
            %             convergence.
            %
            % RETURNS:
            %   conv  - Array of convergence flags for all tests.
            %   vals  - Values checked to obtain `conv`.
            %   names - The names of the convergence checks/equations.
            %   eval  - Logical array into problem.equations indicating which
            %           residual equations we have actually checked 
            %           convergence for.
            %   

            [values, names, tolerances, evaluated] = getConvergenceValuesWells(model, problem);
            convergence = values < tolerances;
        end
        
         function [values, tolerances, names, evalauted] = getFacilityConvergenceValues(model, problem, varargin)
            [values, names, tolerances, evalauted] = getConvergenceValuesWells(model, problem);
         end
  
         function [values, tolerances, names] = getConvergenceValues(model, problem, varargin)
            [values, names, tolerances] = getConvergenceValuesWells(model, problem);
        end


    end

    methods (Static)
        function [wc, varargout] = handleRepeatedPerforatedcells(wc, varargin)
            % Handle multiple wells perforated in the same cells 
            %
            % SYNOPSIS:
            %   [wc, src1, src2] = handleRepeatedPerforatedcells(wc, src1, src2);
            %
            % DESCRIPTION:
            %   This function treats repeated indices in wc (typically due to
            %   multiple wells intersecting a single cell). The output will
            %   have no repeats in wc, and add together any terms in cqs.
            %
            % PARAMETERS:
            %   wc       - Well cells with possible repeats.
            %   varargin - Any number of arrays that should be processed
            %              to account for repeated entries.
            %
            % RETURNS:
            %   wc        - Well cells with repeats removed.
            %   varargout - Variable inputs processed to fix repeated
            %               indices.
            %

            varargout = varargin;
            [c, ic, ic] = uniqueStable(wc);                     %#ok<ASGLU>
            if numel(c) ~= numel(wc)
                A = sparse(ic, (1:numel(wc))', 1, numel(c), numel(wc));
                wc = c;
                for srcNo = 1:numel(varargin)
                    cqs = varargin{srcNo};
                    for k=1:numel(cqs)
                        cqs{k} = A*cqs{k};
                    end
                    varargout{srcNo} = cqs;
                end
            end
        end
        
        function d = combineCellData(data, ix)
            d = cellfun(@(x) x{ix}, data, 'UniformOutput', false);
            d = vertcat(d{:});
        end
    end
end

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

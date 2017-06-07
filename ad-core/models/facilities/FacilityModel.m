classdef FacilityModel < PhysicalModel
    properties
        WellModels

        toleranceWellBHP
        toleranceWellRate
        toleranceWellMS
        ReservoirModel

        VFPTablesInjector
        VFPTablesProducer
    end

    properties (SetAccess = protected)
        % Canonical list of all extra primary variables added by the wells
        addedPrimaryVarNames = {};
        % Canonical list of additional equations
        addedEquationNames = {};
        % Canonical list of the types of the added equations
        addedEquationTypes = {};
    end

    methods
        function model = FacilityModel(reservoirModel, varargin)
            model = model@PhysicalModel([]);

            % Convergence tolerance for BHP-type controls
            model.toleranceWellBHP  = 1*barsa;
            % Convergence tolerance for rate-type controls
            model.toleranceWellRate = 1/day;
            % Convergence tolerance for multisegment wells
            model.toleranceWellMS   = 1e-6;
            % VFP Tables. EXPERIMENTAL.
            model.VFPTablesInjector = {};
            % VFP Tables. EXPERIMENTAL.
            model.VFPTablesProducer = {};
            model = merge_options(model, varargin{:});
            model.ReservoirModel = reservoirModel;
            model.WellModels = {};
        end

        function model = setupWells(model, W, wellmodels)
            % Set up well models for changed controls or the first
            % simulation step.
            %
            % INPUT:
            % 
            % W       - Well struct (obtained from e.g. addWell or
            %           processWells)
            %
            % wellmodels (OPTIONAL ARGUMENT)
            %          - Cell array of equal length to W, containing class
            %          instances for each well (e.g. SimpleWell,
            %          MultisegmentWell, or classes derived from these). 
            %          If not provided, well models be constructed from the
            %          input. 
            nw = numel(W);
            if model.getNumberOfWells == 0
                % First time setup
                [pvars, eqnames, eqtypes] = deal(cell(nw, 1));
                model.WellModels = cell(nw, 1);
                for i = 1:nw
                    % Set up models. SimpleWell for the time being
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
                    wm.dpMaxRel = model.ReservoirModel.dpMaxRel;

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
                    pvars{i} = wm.getExtraPrimaryVariableNames(model.ReservoirModel);
                    [eqnames{i}, eqtypes{i}] = wm.getExtraEquationNames(model.ReservoirModel);
                    model.WellModels{i} = wm;
                end
                % Combine the different equations and types added by the
                % different wells into a canonical ordering.
                model.addedPrimaryVarNames = uniqueStable([pvars{:}]);
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

        function W = getWellStruct(model)
            % Compute number of wells in facility
            W = cellfun(@(x) x.W, model.WellModels, 'UniformOutput', false);
            W = vertcat(W{:});
        end

        function nwell = getNumberOfWells(model)
            % Compute number of wells in facility
            nwell = numel(model.WellModels);
        end
        function actIx = getIndicesOfActiveWells(model)
            % Compute number of wells in facility
            act = model.getWellStatusMask();
            actIx = (1:numel(act))';
            actIx = actIx(act);
        end

        function act = getWellStatusMask(model)
            % Compute number of wells in facility
            act = cellfun(@(x) x.W.status, model.WellModels);
        end

        function names = getPrimaryVariableNames(model)
            % This includes both the basic variables, and the variables
            % added by complex wells (if any)
            names = [model.getBasicPrimaryVariableNames(), model.addedPrimaryVarNames];
        end

        function names = getBasicPrimaryVariableNames(model)
            % Basic primary variables are phase rates + bhp for active
            % phases in the model.
            actPh = model.ReservoirModel.getActivePhases();
            names = {'qWs', 'qOs', 'qGs', 'bhp'};
            names = names([actPh, true]);
        end

        function [variables, names, map] = getBasicPrimaryVariables(model, wellSol)
            % Get phase rates for the active phases and the bhp.
            % In addition, the map contains indicators used to
            % find the phase rates and BHP values in "variables"
            % since these are of special importance to many
            % applications and are considered canonical (i.e. they
            % are always solution variables in MRST, and functions
            % can assume that they will always be found in the
            % variable set for wells).
            
            if model.getNumberOfWells() == 0
                [variables, names] = deal({});
                [isBHP, isRate] = deal([]);
            else
                % Take the active wellSols
                active = model.getWellStatusMask();
                wellSol = wellSol(active);
                % Get number of phases
                actPh = model.ReservoirModel.getActivePhases();
                bhp = vertcat(wellSol.bhp);
                qWs = vertcat(wellSol.qWs);
                qOs = vertcat(wellSol.qOs);
                qGs = vertcat(wellSol.qGs);
                rates = {qWs, qOs, qGs};
                rates = rates(actPh);
                variables = [rates, bhp];
                names = model.getBasicPrimaryVariableNames();
                
                isBHP = false(size(variables));
                isBHP(end) = true;
                isRate = ~isBHP;
            end
            map = struct('isBHP', isBHP, 'isRate', isRate);
        end

        function [variables, names, wellmap] = getExtraPrimaryVariables(model, wellSol)
        % Extra primary variables are variables required by more
        % advanced wells that are in addition to the basic facility
        % variables (rates + bhp).
        %
        % OUTPUT:
        %
        % names     - Column of cells. Each cell is a string with the name
        %             of an extra-variable.
        %
        % variables - Column of cells. Each element, variables{i}, is a vector given
        %             the value corresponding to extra-variable with name
        %             names{i}. This vector is composed of stacked values
        %             over all the wells that contains this extra-variable.
        %
        % wellmap   - The facility model contains the extra-variables of all
        %             the well models that are used. Let us consider the well
        %             with well number wno (in the set of active wells), then
        %             the Well model is belongs to has its own
        %             extra-variables (a subset of those of the Facility
        %             model). We consider the j-th extra-variable of
        %             the Well model. Then, i = extraMap(wno, j)
        %             says that this extra-variable corresponds to
        %             names{i}.

            actIx = model.getIndicesOfActiveWells();
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

            [basic, bnames, wellMap] = model.getBasicPrimaryVariables(wellSol);
            [extra, enames, extraMap] = model.getExtraPrimaryVariables(wellSol);
            
            wellMap.isBHP = [wellMap.isBHP, false(size(extra))];
            wellMap.isRate = [wellMap.isRate, false(size(extra))];
            
            wellMap.extraMap = extraMap;
            names = [bnames, enames];
            variables = [basic, extra];
        end

        function [sources, wellSystem, wellSol] = getWellContributions(model, wellSol0, wellSol, wellvars, wellMap, p, mob, rho, dissolved, comp, dt, iteration)
            % Get the source terms due to the wells, control and well
            % equations and updated well sol. Main gateway for adding wells
            % to a set of equations.
            if isnan(iteration) || iteration < 0
                warning(['Iteration number is not passed on to well model,', ...
                         'this may indicate wellbore pressure-drop will never be updated']);
            end
            actWellIx = model.getIndicesOfActiveWells();
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
                varmaps{varNo} = model.getWellVariableMap(addedVars{varNo});
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
                srcMass{phNo} = combineCellData(allMass, phNo);
                srcVol{phNo} = combineCellData(allVol, phNo);
                eqs{phNo} = combineCellData(allBaseEqs, phNo);
            end
            % Components are ordered canonically by reservoir model
            srcComp = cell(1, ncomp);
            for cNo = 1:ncomp
                srcComp{cNo} = combineCellData(allComp, cNo);
            end
            % If we have extra equations, add them in
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
            ctrleq = vertcat(allCtrl{:});

            wc = model.getActiveWellCells();
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

        function wellSol = updateWellSolAfterStep(model, wellSol)
            % Figure out if wells are shut, or changed ontrols
            for wno = 1:numel(wellSol)
                wm = model.WellModels{wno};
                wellSol(wno) = wm.updateWellSolAfterStep(model.ReservoirModel, wellSol(wno));
            end
        end

        function wc = getWellCells(model)
            c = cellfun(@(x) x.W.cells, model.WellModels, 'UniformOutput', false);
            wc = vertcat(c{:});
        end

        function wc = getActiveWellCells(model)
            c = cellfun(@(x) x.W.cells, model.WellModels, 'UniformOutput', false);
            active = model.getWellStatusMask();
            wc = vertcat(c{active});
        end

        function ws = setWellSolStatistics(model, ws, sources)
            % Store extra output, typically black oil-like
            p2w = getPerforationToWellMapping(model.getWellStruct());
            % Map into active wells
            active = model.getWellStatusMask();
            p2w = p2w(active(p2w));

            gind = model.ReservoirModel.getPhaseIndex('G');
            oind = model.ReservoirModel.getPhaseIndex('O');
            wind = model.ReservoirModel.getPhaseIndex('W');
            srcRes = cellfun(@double, sources.phaseVolume, 'UniformOutput', false);
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

                ws(i).flux = sum(qR(p2w == i, :), 2);
            end
        end

        % Implementation details for stand-alone model
        function [wellSol, restVars] = updateWellSol(model, wellSol, problem, dx, drivingForces, restVars) %#ok
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
            actIx = model.getIndicesOfActiveWells();
            nW = numel(actIx);
            activeVars = false(nW, nVars);
            dx_well = cell(nW, nVars);
            % Loop over all possible well variables
            for varNo = 1:nVars
                wf = wellVars{varNo};
                if any(strcmp(restVars, wf))
                    dv = model.getIncrement(dx, problem, wf);

                    isVarWell = model.getWellVariableMap(wf);
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

        function isVarWell = getWellVariableMap(model, wf)
            act = model.getWellStatusMask();
            counts = cellfun(@(x) x.getVariableCounts(wf), model.WellModels(act));
            isVarWell = rldecode((1:nnz(act))', counts);
        end

        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
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
            if isfield(state, 'wellSol')
                state.wellSol = model.updateWellSol(state.wellSol, problem, dx, drivingForces);
                % Handle the directly assigned values (i.e. can be deduced directly from
                % the well controls.
                W = drivingForces.W;
                phIndices = model.ReservoirModel.getPhaseIndices();
                state.wellSol = assignWellValuesFromControl(model.ReservoirModel, state.wellSol, W, phIndices(1), phIndices(2), phIndices(3));
            end
            report = [];
        end

        function state = validateState(model, state)
            if ~isfield(state, 'wellSol') || isempty(state.wellSol),
                if isfield(state, 'wellSol'),
                    state = rmfield(state, 'wellSol');
                end
                W = model.getWellStruct();
                state.wellSol = initWellSolAD(W, model.ReservoirModel, state);
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
            forces = getValidDrivingForces@PhysicalModel(model);
            forces.W   = [];
            forces.bc  = [];
            forces.src = [];
        end

        function [convergence, values, names, evaluated] = checkFacilityConvergence(model, problem)
            % For checking on the subset of variables specific to the
            % facility
            [convergence, values, evaluated, names] = checkWellConvergence(model, ...
                                                              problem);

        end

        function [convergence, values, names] = checkConvergence(model, problem, varargin)
            % Used when facility is run as a stand-alone model
            [convergence, values, names] = ...
                model.checkFacilityConvergence(problem);
        end
    end

    methods (Static)
        function [wc, varargout] = handleRepeatedPerforatedcells(wc, varargin)
            % This function treats repeated indices in wc (typically due to
            % multiple wells intersecting a single cell). The output will
            % have no repeats in wc, and add together any terms in cqs.
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
    end
end

function d = combineCellData(data, ix)
    d = cellfun(@(x) x{ix}, data, 'UniformOutput', false);
    isAD = ~cellfun(@isnumeric, d);
    if any(~isAD) && any(isAD)
        s = d(isAD);
        s = s{1};
        for i = 1:numel(isAD)
            if ~isAD(i)
                d{i} = double2ADI(d{i}, s);
            end
        end
    end
    d = vertcat(d{:});
end

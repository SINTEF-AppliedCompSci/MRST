classdef FacilityModel < PhysicalModel
    properties
        WellModels
        
        toleranceWellBHP
        toleranceWellRate
        ReservoirModel
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
            
            model.toleranceWellBHP = 1*barsa;
            model.toleranceWellRate = 1/day;
            
            model = merge_options(model, varargin{:});
            model.ReservoirModel = reservoirModel;
            model.WellModels = {};
        end
        
        function model = setupWells(model, W, wellmodels)
            % Set up wells for changed controls or first simulation step
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
                    model.WellModels{i}.updateWell(W(i));
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

        function [rates, bhp, names] = getBasicPrimaryVariables(model, wellSol)
            % Get phase rates + bhp for active phases
            actPh = model.ReservoirModel.getActivePhases();
            bhp = vertcat(wellSol.bhp);
            qWs = vertcat(wellSol.qWs);
            qOs = vertcat(wellSol.qOs);
            qGs = vertcat(wellSol.qGs);
            rates = {qWs, qOs, qGs};
            rates = rates(actPh);
            
            names = model.getBasicPrimaryVariableNames();
        end
        
        function [variables, names, wellmap] = getExtraPrimaryVariables(model, wellSol)
            % Extra primary variables are variables required by more
            % advanced wells that are in addition to the basic facility
            % variables (rates + bhp).
            names = model.addedPrimaryVarNames;
            nw = model.getNumberOfWells();
            nv = numel(names);
            vars = cell(nw, nv);
            
            wellmap = zeros(nw, nv);
            if nv > 0
                all_ix = (1:nv)';
                for i = 1:nw
                    [v, n] = model.WellModels{i}.getExtraPrimaryVariables(wellSol(i), model.ReservoirModel);

                    for j = 1:numel(v)
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
        
        function [srcMass, srcVol, eqs, ctrleq, names, types, wellSol] = getWellContributions(model, wellSol0, wellSol, qWell, bhp, wellvars, wellMap, p, mob, rho, comp, dt, iteration)
            % Get the source terms due to the wells, control and well
            % equations and updated well sol. Main gateway for adding wells
            % to a set of equations.
            if isnan(iteration) || iteration < 0
                warning(['Iteration number is not passed on to well model,', ...
                         'this may indicate wellbore pressure-drop will never be updated']);
            end
            
            nw = model.getNumberOfWells();
            
            allBaseEqs = cell(nw, 1);
            allCtrl = cell(nw, 1);
            allVol = cell(nw, 1);
            allMass = cell(nw, 1);
            
            enames = model.addedEquationNames;
            etypes = model.addedEquationTypes;
            n_extra = numel(enames);
            assert(numel(etypes) == n_extra);
            
            allExtraEqs = cell(nw, n_extra);

            
            addedVars = model.addedPrimaryVarNames;
            maps = cell(1, numel(addedVars));
            for varNo = 1:numel(addedVars)
                maps{varNo} = model.getWellVariableMap(addedVars{varNo});
            end
            
            [basenames, basetypes] = model.WellModels{1}.getWellEquationNames(model.ReservoirModel);
            for i = 1:nw
                wm = model.WellModels{i};

                W = wm.W;
                wc = W.cells;
                pw = p(wc);
                mobw = getCellSubset(mob, wc);
                rhow = getCellSubset(rho, wc);
                compw = getComponentCellSubset(comp, wc);
                varw = getVariableSubsetWell(wellvars, maps, i);
                
                % Renumber to the ordering of variables for this well
                renum = wellMap(i, wellMap(i, :) > 0);
                varw = varw(renum);
                
                qw = cellfun(@(x) x(i), qWell, 'uniformoutput', false);
                bh = bhp(i);
                % Update pressure
                wellSol(i) = wm.updateConnectionPressureDrop(wellSol0(i), wellSol(i), model.ReservoirModel, qw, bh, varw, pw, mobw, rhow, compw, dt, iteration);
                % Update limits
                [qw, bh, wellSol(i), ok] = wm.updateLimits(wellSol0(i), wellSol(i), model.ReservoirModel, qw, bh, varw, pw, mobw, rhow, compw, dt, iteration);
                if ~ok
                    bhp(i) = bh;
                    for phNo = 1:numel(qw)
                        qWell{phNo}(i) = qw{phNo};
                    end
                end
               % Set up well equations and source terms
               [allBaseEqs{i}, allCtrl{i}, extraEqs, extraNames, allMass{i}, allVol{i}, wellSol(i)] =...
                   wm.computeWellEquations(wellSol0(i), wellSol(i), model.ReservoirModel, qw, bh, varw, pw, mobw, rhow, compw, dt, iteration);
               for eqNo = 1:numel(extraEqs)
                   % Map into global list of equations
                   ix = strcmpi(enames, extraNames{eqNo});
                   allExtraEqs{i, ix} = extraEqs{eqNo};
               end
            end
            % We have assembled all equations for each well. Combine the
            % equations from the different wells into one (array) of each
            % type.
            nPh = nnz(model.ReservoirModel.getActivePhases);
            [srcMass, srcVol, eqs] = deal(cell(1, nPh));
            for phNo = 1:nPh
                srcMass{phNo} = combineCellData(allMass, phNo);
                srcVol{phNo} = combineCellData(allVol, phNo);
                eqs{phNo} = combineCellData(allBaseEqs, phNo);
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
        end

        function wellSol = updateWellSolAfterStep(model, resmodel, wellSol)
            % Figure out if wells are shut, or changed ontrols
            for wno = 1:numel(wellSol)
                wm = model.WellModels{wno};
                wellSol(wno) = wm.updateWellSolAfterStep(resmodel, wellSol(wno));
            end
        end
        
        function wc = getWellCells(model)
            c = cellfun(@(x) x.W.cells, model.WellModels, 'UniformOutput', false);
            wc = vertcat(c{:});
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
            resModel = model.ReservoirModel;
            for i = 1:numel(wellVars)
                wf = wellVars{i};
                dv = model.getIncrement(dx, problem, wf);

                if strcmpi(wf, 'bhp')
                    % Bottom hole is a bit special - we apply the pressure update
                    % limits here as well.
                    bhp = vertcat(wellSol.bhp);
                    dv = resModel.limitUpdateRelative(dv, bhp, resModel.dpMaxRel);
                    dv = resModel.limitUpdateAbsolute(dv, resModel.dpMaxAbs);
                end
                isVarWell = model.getWellVariableMap(wf);
                for j = 1:numel(wellSol)
                    subs = isVarWell == j;
                    if any(subs)
                        % Put this into WellModel so that we can have
                        % limits on properties
                        wellSol(j) = model.WellModels{j}.incrementProp(wellSol(j), wf, dv(subs));
                    end
                end
                % Field is taken care of
                restVars = model.stripVars(restVars, wf);
            end
        end
        
        function isVarWell = getWellVariableMap(model, wf)
            nw = model.getNumberOfWells();
            counts = cellfun(@(x) x.getVariableCounts(wf), model.WellModels);
            isVarWell = rldecode((1:nw)', counts);
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            opt = struct('iteration', nan, 'resOnly', false);
            opt = merge_options(opt, varargin{:});
            wellSol = state.wellSol;
            wellSol0 = state0.wellSol;
            resmodel = model.ReservoirModel;
            % Get variables from facility and wells
            [qWell, bhp, basicWellNames] = model.getBasicPrimaryVariables(wellSol);
            [wellVars, wellExtraNames, wellMap] = model.getExtraPrimaryVariables(wellSol);
            if ~opt.resOnly
                [qWell{:}, bhp, wellVars{:}] = initVariablesADI(qWell{:}, bhp, wellVars{:});
            end
            
            if isa(resmodel, 'ThreePhaseBlackOilModel')
                [rs, rv] = resmodel.getProps(state, 'rs', 'rv');
            else
                [rs, rv] = deal([]);
            end
            assert(isfield(state, 'rho'));
            assert(isfield(state, 'mob'));
            p = resmodel.getProp(state, 'pressure');
            
            nPh = size(state.rho, 2);
            [mob, rho] = deal(cell(1, nPh));
            for i = 1:nPh
                mob{i} = state.mob(:, i);
                rho{i} = state.rho(:, i);
            end
            components = resmodel.getDissolutionMatrix(rs, rv);
            [srcMass, srcVol, weqs, ctrleq, wnames, wtypes, state.wellSol] = ...
                model.getWellContributions(wellSol0, wellSol, qWell, bhp, wellVars, ...
                        wellMap, p, mob, rho, components, dt, opt.iteration);
            
            eqs = {weqs{:}, ctrleq};
            names = {wnames{:}, 'closureWells'};
            types = {wtypes{:}, 'well'};
            
            primaryVars = {basicWellNames{:}, wellExtraNames{:}};
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            if isfield(state, 'wellSol')
                state.wellSol = model.updateWellSol(state.wellSol, problem, dx, drivingForces);
            end
            report = [];
        end
        
        function state = validateState(model, state)
            if ~isfield(state, 'wellSol')
                if ~isfield(state, 'wellSol') || isempty(state.wellSol),
                   if isfield(state, 'wellSol'),
                      state = rmfield(state, 'wellSol');
                   end
                   W = model.getWellStruct();
                   state.wellSol = initWellSolAD(W, model.ReservoirModel, state);
                end
            end
            for wno = 1:numel(model.WellModels)
                new_ws = model.WellModels{wno}.validateWellSol(model.ReservoirModel, state.wellSol(wno));
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
            [convergence, values, evaluated, names] = checkWellConvergence(model, problem);
        end
        
        function [convergence, values, names] = checkConvergence(model, problem, varargin)
            % Used when facility is run as a stand-alone model
            [convergence, values, names] = model.checkFacilityConvergence(problem);
        end
    end
end

function celldata = getComponentCellSubset(celldata, wc)
    for i = 1:numel(celldata)
        for j = 1:numel(celldata{i});
            if ~isempty(celldata{i}{j})
                celldata{i}{j} = celldata{i}{j}(wc);
            end
        end
    end
end

function d = combineCellData(data, ix)
    d = cellfun(@(x) x{ix}, data, 'UniformOutput', false);
    d = vertcat(d{:});
end

function subset = getCellSubset(celldata, wc)
    subset = cell(size(celldata));
    for i = 1:numel(subset)
        if ~isempty(celldata{i})
            subset{i} = celldata{i}(wc);
        end
    end
end

function subset = getVariableSubsetWell(vars, wellMap, ix)
    subset = cell(size(vars));
    for i = 1:numel(subset)
        subset{i} = vars{i}(wellMap{i} == ix);
    end
end
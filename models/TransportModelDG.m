classdef TransportModelDG < TransportModel
    
    properties
        discretization = []    % Discretization
        dgVariables    = {'s'} % Transport variables we discretize with dG
        limiters               % Limiters
        storeUnlimited = false % Store unlimited state for plotting/debug
        dsMaxAbsDiv    = 4     % Control reduction of parentModel.dsMaxAbs
    end
    
    methods
        %-----------------------------------------------------------------%
        function model = TransportModelDG(parent, varargin)
            % Parent model initialization
            model = model@TransportModel(parent); 
            % Default limiters
            names    = model.dgVariables;
            limits   = {[0,1]}; % Limiter limits
            tol      = 0;       % Limiter tolerances
            limiters = [];      % Add limiters
            limiters = addLimiter(limiters           , ... % TVB limiter
                                  'type'     , 'tvb' , ...
                                  'variables', names , ...
                                  'limits'   , limits, ...
                                  'tol'      , tol   );
            limiters = addLimiter(limiters            , ... % Scale limiter
                                  'type'     , 'scale', ...
                                  'variables', names  , ...
                                  'limits'   , limits , ...
                                  'tol'      , tol    );
            model.limiters       = limiters;
            model.storeUnlimited = false;
            % Merge options
            [model, discretizationArgs] = merge_options(model, varargin{:});
            % Construct discretization
            if isempty(model.discretization)
                model.discretization = DGDiscretization(model.G, discretizationArgs{:});
            end
            % Get limiters
            for l = 1:numel(model.limiters)
                limiter = model.limiters(l);
                model.limiters(l).function = getLimiter(model, limiter.type);
            end
            % Set up DG operators
            model.parentModel.operators = setupOperatorsDG(model.discretization  , ...
                                                           model.parentModel.G   , ...
                                                           model.parentModel.rock);
            % Assign discretization to parentModel (currently stored on
            % operators)
            model.parentModel.operators.discretization = model.discretization;
            % Pressure is not solved with DG, make sure to don't store
            % things to state that should be recomputed in pressure step
            model.parentModel.outputFluxes         = false;
            model.parentModel.OutputStateFunctions = {};
            model.parentModel.useCNVConvergence  = false;
            model.parentModel.nonlinearTolerance = 1e-3;
        end
        
        %-----------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name, varargin)
            % Get variable fiels, check if it is dof
            isDof = any(strcmpi(name(1:end-3), model.dgVariables));
            if isDof
                lookup = name(1:end-3);
            else
                lookup = name;
            end
            [fn, index] = getVariableField@TransportModel(model, lookup, varargin{:});
            if isDof && ~isempty(fn)
                fn = [fn, 'dof'];
            end
        end
        
        %-----------------------------------------------------------------%
        function groupings = getStateFunctionGroupings(model)
            groupings = model.parentModel.getStateFunctionGroupings();
        end
        
        %-----------------------------------------------------------------%
        function state = validateState(model, state)
            % Set degree in each cell
            state.degree = repmat(model.discretization.degree, model.G.cells.num, 1);
            % Well are treated as dG(0)
            wm = model.parentModel.FacilityModel.WellModels;
            for i = 1:numel(wm)
                state.degree(wm{i}.W.cells,:) = 0;
            end
            % Let parent model do its thing
            state = validateState@TransportModel(model, state);
            % Assign dofs
            state = assignDofFromState(model.discretization, state, model.dgVariables);
            % Evaluate basis functions in cubature points
            state = model.discretization.evaluateBasisFunctions(state, inf, inf);
        end
        
        %-----------------------------------------------------------------%
        function [vars, names, origin] = getPrimaryVariables(model, state)
            % Get primary variables
            [vars, names, origin] = model.parentModel.getPrimaryVariables(state);
            isParent = strcmp(origin, class(model.parentModel));
            vars     = vars(isParent);
            names    = names(isParent);
            % If saturation is a dG variable, we have already assigned relace 's' by
            % 'sW', etc
            isSat = strcmpi(model.dgVariables, 's');
            model.dgVariables(isSat) = [];
            % Add dof to ending of dG variable names
            isDof        = ismember(names, model.dgVariables);
            names(isDof) = cellfun(@(bn) [bn, 'dof'], names(isDof), 'UniformOutput', false);
            % Replace variables with dofs
            origin = origin(isParent);
            isBO   = strcmpi(origin, 'GenericBlackOilModel');
            for i = 1:numel(names)
                if isDof(i)
                    [fn, ~] = model.getVariableField(names{i}, false);
                    if ~isempty(fn)
                        vars{i} = model.getProp(state, names{i});
                    elseif any(strcmpi(names{i}, {'x', 'xdof'})) && isBO(i)
                        error('dG currently does not support disgas/vapoil')
                    end
                end
            end
        end
        %-----------------------------------------------------------------%
        function [state, names, origin] = getStateAD(model, state, init)
            % Get AD state (handled by getStateDG)
            [cells, faces] = deal(Inf);
            [state, names, origin] = model.getStateDG(state, cells, faces, init);
        end
        %-----------------------------------------------------------------%
        function [state, names, origin] = getStateDG(model, state, cells, faces, init)
            if nargin < 5
                init = true;
            end            
            parent = model.parentModel;
            % Get the AD state for this model
            [basevars, basenames, baseorigin] = model.getPrimaryVariables(state);
            isParent   = strcmp(baseorigin, class(parent));
            basevars   = basevars(isParent);
            basenames  = basenames(isParent);
            baseorigin = baseorigin(isParent);
            % Find saturations
            isS = false(size(basevars));
            nph = parent.getNumberOfPhases();
            phase_variable_index = zeros(nph, 1);
            for i = 1:numel(basevars)
                [f, ix] = model.getVariableField(basenames{i}, false);
                if any(strcmpi(f, {'s', 'sdof'}))
                    isS(i) = true;
                    phase_variable_index(ix) = i;
                end
            end
            % Figure out saturation logic
            isP    = strcmp(basenames, 'pressure');
            vars   = basevars;
            names  = basenames;
            origin = baseorigin;
            useTotalSaturation = ....
                strcmpi(model.formulation, 'totalSaturation') && sum(isS) == nph - 1;
            if useTotalSaturation
                % Replace pressure with total saturation
                if any(strcmpi('sT', model.dgVariables))
                    replacement = 'sTdof';
                else
                    replacement = 'sT';
                end
                sTdof       = model.getProp(state, replacement);
                % Replacing
                vars{isP}   = sTdof;
                names{isP}  = replacement;
                origin{isP} = class(model);
            else
                % Remove pressure and skip saturation closure
                vars   = vars(~isP);
                names  = names(~isP);
                origin = origin(~isP);
            end
            if init
                [vars{:}] = model.AutoDiffBackend.initVariablesAD(vars{:});
            end
            if useTotalSaturation
                basevars(~isP) = vars(~isP);
            else
                basevars(~isP) = vars;
            end
            % Evluate basis functions for use later
            state = model.discretization.evaluateBasisFunctions(state, cells, faces);
            isDof = false(size(names));
            [cellMean, cellVars, faceVars] = deal(cell(size(vars)));
            for i = 1:numel(vars)
                [cellMean{i}, cellVars{i}, faceVars{i}, isDof(i)] ...
                    = model.evaluateBaseVariable(state, basevars{i}, basenames{i});
            end
            % Let parent model handle initStateAD
            basenames(isDof) = cellfun(@(bn) bn(1:end-3), basenames(isDof), 'UniformOutput', false);
            % Initialize cell mean state
            if ~all(isinf(state.mcells))
                parent.G.cells.num = numel(state.mcells);
            end
            state = parent.initStateAD(state, cellMean, basenames, baseorigin);
            % First, store dofs to state (only used in case we have BCs)
            state = model.assignDofsToState(state, vars, names);
            % Evaluate non-dg and non-ad variables in cubature points
            state = model.evaluateBaseVariables(state);
            % Initialize well state
            state.wellStateDG = parent.initStateAD(state.wellStateDG, cellMean, basenames, baseorigin);
            % Initialize cell state
            parent.G.cells.num = numel(value(cellVars{1}));
            state.cellStateDG = parent.initStateAD(state.cellStateDG, cellVars, basenames, baseorigin);
            % Initialize face state
            parent.G.cells.num = numel(value(faceVars{1}));
            state.faceStateDG = parent.initStateAD(state.faceStateDG, faceVars, basenames, baseorigin);
            if useTotalSaturation
                % Set total saturation as well
                sTdof       = vars{isP};
                state.sTdof = sTdof;
                [meanValue, cellValue, faceValue] = model.evaluateBaseVariable(state, sTdof, 'sTdof');
                state.wellStateDG = model.setProp(state.wellStateDG, 'sT', meanValue);
                state.cellStateDG = model.setProp(state.cellStateDG, 'sT', cellValue);
                state.faceStateDG = model.setProp(state.faceStateDG, 'sT', faceValue);
                state             = model.setProp(state, 'sT', meanValue);
            end
        end
        
        %-----------------------------------------------------------------%
        function state = assignDofsToState(model, state, vars, names)
            % Assign dofs to state for e.g., BC evaluation
            fill = zeros(sum(state.nDof), 1);
            ix   = model.discretization.getDofIx(state, 1);
            fill(ix) = 1;
            snames = model.parentModel.getSaturationVarNames();
            snames = cellfun(@(sn) [sn, 'dof'], snames, 'UniformOutput', false);
            sdof   = cell(1, numel(snames));
            removed  = false(size(vars));
            for i = 1:numel(vars)
                ix = strcmpi(names{i}, snames);
                if any(ix)
                    sdof{ix} = vars{i};
                    fill = fill - vars{i};
                    removed(i) = true;
                end
            end
            ix = strcmpi(setdiff(snames, lower(names)), snames);
            sdof{ix} = fill;
            state = model.setProp(state, 'sdof', sdof);
            
            vars  = vars(~removed);
            names = names(~removed);
            for i = 1:numel(vars)
                state = model.setProp(state, names{i}, vars{i});
            end
            
        end
        
        %-----------------------------------------------------------------%
        function [meanVal, cellVal, faceVal, isDof] = evaluateBaseVariable(model, state, var, name)
            assert(isfield(state, 'psi_c') && isfield(state, 'psi_c'));
            isDof = false;
            if any(strcmpi(name(1:end-3), model.dgVariables)) && strcmpi(name(end-2:end), 'dof')
                % dG - do evaluation at cubature points
                isDof   = true;
                meanVal = model.discretization.getCellMean(state, state.mcells, var);
                cellVal = model.discretization.evaluateProp(state, var, 'cell', state.mcells);
                faceVal = model.discretization.evaluateProp(state, var, 'face', state.mfaces);
            else
                % Not dG - repeat to match number of cubature points
                c = state.mcells;
                if isinf(state.mcells), c = ':'; end
                meanVal = var(c,:);
                cellVal = var(state.cells,:);
                faceVal = var(state.fcells,:);
            end
        end
        
        %-----------------------------------------------------------------%
        function state = evaluateBaseVariables(model, state)
            % Evaluate all non-dG and non-primary variables at all cubature
            % points
            [cellStateDG, faceStateDG, wellStateDG] = deal(state);
            if ~(isfield(state, 'psi_c') && isfield(state, 'psi_f'))
                % Evaluate basis functions at cubature points
                state = model.discretization.evaluateBasisFunctions(state, Inf, Inf);
            end
            % Assign type and cells/faces
            cellStateDG.type   = 'cell';
            cellStateDG.cells  = state.cells;
            cellStateDG.fcells = state.fcells;
            cellStateDG.faces  = state.faces;
            wellStateDG.type   = 'cell';
            wellStateDG.cells  = state.cells;
            faceStateDG.type   = 'face';
            faceStateDG.cells  = state.fcells;
            faceStateDG.faces  = state.faces;
            % Evaluate valriables
            names = fieldnames(state);
            for k = 1:numel(names)
                name = names{k};
                [fn , index] = model.getVariableField(name, false);
                if ~isempty(fn) && isa(state.(fn), 'double') && ... % Only doubles
                        ~any(strcmpi(name, model.dgVariables))        % ... dG variables are set from dofs
                    % ... and only variables of correct dimension
                    n = size(double(state.(fn)),1);
                    if (n ~= model.G.cells.num && n ~= sum(state.nDof))
                        continue
                    else
                        % Evaluate
                        [meanVal, cellVal, faceVal] ...
                            = model.evaluateBaseVariable(state, state.(name)(:,index), name);
                        if strcmpi(name(end-2:end), 'dof')
                            name = name(1:end-3);
                        end
                        % Assign to state
                        cellStateDG = model.setProp(cellStateDG, name, cellVal);
                        wellStateDG = model.setProp(wellStateDG, name, meanVal);
                        faceStateDG = model.setProp(faceStateDG, name, faceVal);
                    end
                end
            end
            % Set flag (compositional models)
            if isfield(faceStateDG, 'flag')
                faceStateDG.flag = faceStateDG.flag(fcells);
            end
            % Store cell/well/face states to state
            state.cellStateDG = cellStateDG;
            state.wellStateDG = wellStateDG;
            state.faceStateDG = faceStateDG;
        end
        
        %-----------------------------------------------------------------%
        function state = assignCellMean(model, state)
            % Assign cell mean for all dg variables
            names = model.dgVariables;
            for name = names
                if isfield(state, [name{1}, 'dof'])
                    dof = model.getProp(state, [name{1}, 'dof']);
                    v   = model.discretization.getCellMean(state, Inf, dof);
                    state.(name{1}) = v;
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function model = validateModel(model, varargin)
            % Validate model before simulation
            isSat = strcmpi('s', model.dgVariables);
            assert(any(isSat), 'dG currently requires that saturation is a dG variable');
            % Add phase saturations as dgVariables
            if isSat
                phNames = model.parentModel.getPhaseNames();
                for ph = phNames
                    model.dgVariables{end+1} = ['s', ph];
                end
                if strcmpi(model.formulation, 'totalSaturation')
                    model.dgVariables{end+1} = 'sT';
                end
            end
            % Parent class handles most of the validation
            model = validateModel@TransportModel(model, varargin{:});
        end
        
        %-----------------------------------------------------------------%
        function model = setupStateFunctionGroupings(model, default)
            if ~default
                % State functions have already been set up
                assert(isa(model.parentModel.FlowDiscretization, 'FlowDiscretizationDG'));
                return
            end
            model = setupStateFunctionGroupings@TransportModel(model, default);
            pModel = model.parentModel;
            % Adjust flow property functions
            fp = model.parentModel.FlowPropertyFunctions;
            fp = fp.setStateFunction('GravityPermeabilityGradient', GravityPermeabilityGradientDG(pModel));
            fp = fp.setStateFunction('TotalMobility', TotalMobility(pModel));
            model.parentModel.FlowPropertyFunctions = fp;
            % Adjust PVT property functions
            pvt    = model.parentModel.PVTPropertyFunctions;
            pvtReg = pvt.getRegionPVT(pModel);
            if isa(pvt.PoreVolume, 'PoreVolume')
                pvt.PoreVolume = PoreVolumeDG(pModel, pvtReg);
            elseif isa(pvt.PoreVolume, 'BlackOilPoreVolume')
                pvt.PoreVolume = BlackOilPoreVolumeDG(pModel, pvtReg);
            else
                error('Pore volume type not supported by dG')
            end
            model.parentModel.PVTPropertyFunctions = pvt;
            % Adjust flux discretization
            model.parentModel.FlowDiscretization = FlowDiscretizationDG(pModel);
        end
    
        %-----------------------------------------------------------------%
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces)
            % Get model equations
            state0 = model.evaluateBaseVariables(state0);
            pmodel = model.parentModel;
            [acc, flux, cellflux, names, types] = pmodel.FlowDiscretization.componentConservationEquations(pmodel, state, state0, dt);
            state.wellStateDG = rmfield(state.wellStateDG, 'FlowProps');
            state.wellStateDG = rmfield(state.wellStateDG, 'FluxDisc');
            src = pmodel.FacilityModel.getComponentSources(state.wellStateDG);
            % Treat source or bc terms
            if ~isempty(drivingForces.bc) || ~isempty(drivingForces.src)
                fluxBC  = model.computeBoundaryConditions(state, state0, dt, drivingForces.bc);
            end
            % Assemble equations and add in sources
            if strcmpi(model.formulation, 'missingPhase')
                % Skip the last phase! Only mass-conservative for
                % incompressible problems
                acc   = acc(1:end-1);
                flux  = flux(1:end-1);
                names = names(1:end-1);
                types = types(1:end-1);
            end
            d        = model.discretization;
            d.nDof   = state.nDof;
            d.dofPos = state.dofPos;
            psi      = d.basis.psi;
            gradPsi  = d.basis.gradPsi;
            ixw      = d.getDofIx(state, 1, src.cells);
            ix       = d.getDofIx(state, Inf);
            d.sample = state.cellStateDG.s{1}(ix)*0;
            eqs      = cell(1, numel(acc));
            state.wellStateDG.cells = (1:pmodel.G.cells.num)';
            
            cells  = rldecode((1:pmodel.G.cells.num)', state.nDof, 1);
            pv     = pmodel.operators.pv(cells);
            rhoS   = pmodel.getSurfaceDensities();
            cnames = pmodel.getComponentNames();
            
            for i = 1:numel(acc)
                eqs{i} = d.inner(acc{i}     , psi    , 'dV') ...
                       - d.inner(cellflux{i}, gradPsi, 'dV') ...
                       + d.inner(flux{i}    , psi    , 'dS');
                if ~isempty(drivingForces.bc)
                    eqs{i} = eqs{i} + d.inner(fluxBC{i}, psi, 'dS', drivingForces.bc.face);
                end
                if ~isempty(src.cells)
                    eqs{i}(ixw) = eqs{i}(ixw) - src.value{i};
                end
                if ~pmodel.useCNVConvergence
                    sub = strcmpi(names{i}, cnames);
                    eqs{i} = eqs{i}.*(dt./(pv.*rhoS(sub)));
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function q = computeBoundaryConditions(model, state, state0, dt, bc)
            % Compute values at boundary faces with prescribed BCs
            % Set cells/faces where we need to compute properties
            faces = bc.face;
            cells = sum(model.G.faces.neighbors(faces,:),2);
            % Compute preoperites with AD
            bcStateDG  = model.getStateDG(state , cells, faces, false);
            % Compute properites without AD (using getStateDG ensures
            % that properties are evaluated consistently)
%             bcStateDG0 = model.getStateDG(state0, cells, faces, false);
            % Get flow state (either bcStateDG or bcStateDG0)
            fsb = FlowStateBuilderDG();
            bcStateDG = fsb.build(model.parentModel.FlowDiscretization, model.parentModel, bcStateDG, [], []);
%             bcStateDG  = model.parentModel.FluxDiscretization.buildFlowState(model.parentModel, bcStateDG, bcStateDG0, dt);
            % Compute boundary condition fluxes
            model.parentModel = model.parentModel.FlowDiscretization.expandRegions(model.parentModel, 'cells', bcStateDG.cells);
            q = computeBoundaryFluxesDG(model.parentModel, bcStateDG, bc);
        end
        
        %-----------------------------------------------------------------%
        function [model, state] = prepareTimestep(model, state, state0, dt, drivingForces)
            % TransportModel handles prepareTimestep
            [model, state] = prepareTimestep@TransportModel(model, state, state0, dt, drivingForces);
        end
        
        %-----------------------------------------------------------------%
        function [restVars, satVars, wellVars] = splitPrimaryVariables(model, vars)
            % Split variables
            vars = cellfun(@(n) n(1:end-3), vars, 'UniformOutput', false);
            [restVars, satVars, wellVars] = model.parentModel.splitPrimaryVariables(vars);
            restVars = cellfun(@(n) [n, 'dof'], restVars, 'UniformOutput', false);
            satVars  = cellfun(@(n) [n, 'dof'], satVars , 'UniformOutput', false);
        end
        
        %-----------------------------------------------------------------%
        function [state, report] = updateState(model, state, problem, dx, drivingForces)  
            % Remove DG states
            state = rmfield(state, 'cellStateDG');
            state = rmfield(state, 'faceStateDG');
            state = rmfield(state, 'wellStateDG');
            % Store state before update
            state0 = state;
            [restVars, satVars] = model.splitPrimaryVariables(problem.primaryVariables);
            % Update saturation dofs
            state = model.updateSaturations(state, dx, problem, satVars);
            % Update non-saturation dofs
            state = model.updateDofs(state, dx, problem, restVars);
            % Update cell averages from dofs
            state = model.assignCellMean(state);
            % Compute dx for cell averages
            dx0 = model.getMeanIncrement(state, state0, problem.primaryVariables);
            % Let parent model do its thing
            problem0 = problem;
            problem0.primaryVariables = cellfun(@(n) n(1:end-3), problem0.primaryVariables, 'UniformOutput', false);
            [state0_corr, report] = updateState@TransportModel(model, state0, problem0, dx0, drivingForces);
            % Correct updates in dofs according to parent model
            dx0_corr = model.getMeanIncrement(state0_corr, state0, problem.primaryVariables);
            cells    = rldecode((1:model.G.cells.num)', state.nDof, 1);
            frac     = cellfun(@(x,y) x(cells)./y(cells), dx0_corr, dx0, 'UniformOutput', false);
            for i = 1:numel(frac)
                frac{i}(~isfinite((frac{i}))) = 1;
            end
            dx_corr  = cellfun(@(dx, f) dx.*f, dx, frac, 'UniformOutput', false);
            % Update saturation dofs
            state = model.updateSaturations(state0, dx_corr, problem, satVars);
            % Update non-saturation dofs
            state = model.updateDofs(state, dx_corr, problem, restVars);
            % Update cell averages from dofs
            state = model.assignCellMean(state);
        end
        
        %-----------------------------------------------------------------%
        function dx = getMeanIncrement(model, state, state0, vars)
            % Get the increment in the mean value of a set of variables
            dx = cell(numel(vars),1);
            for i = 1:numel(vars)
                vn    = vars{i}(1:end-3);
                v     = model.getProp(state, vn);
                v0    = model.getProp(state0, vn);
                dx{i} = v - v0;
            end
        end
        
        % ----------------------------------------------------------------%
        function state = updateDofs(model, state, dx, problem, dofVars)
            % Update non-saturation dofs
            for i = 1:numel(dofVars)
                dvMaxAbs = inf;
                if strcmpi(dofVars{i}, 'sTdof')
                    dvMaxAbs = model.dsMaxTotal./min(model.discretization.basis.nDof,model.dsMaxAbsDiv);
                end
                state = updateStateFromIncrement(model, state, dx, problem, dofVars{i}, inf, dvMaxAbs);
            end
        end
        
        % ----------------------------------------------------------------%
        function state = updateSaturations(model, state, dx, problem, satVars)
            % Update dG saturation variables
            if nargin < 5
                % Get the saturation names directly from the problem
                [~, satVars] = ...
                    splitPrimaryVariables(model, problem.primaryVariables);
            end
            if isempty(satVars)
                % No saturations passed, nothing to do here.
                return
            end
            % Solution variables should be saturations directly, find the
            % missing link
            saturations0 = lower(model.parentModel.getSaturationVarNames);
            saturations  = cellfun(@(n) [n, 'dof'], saturations0, 'uniformOutput', false);
            fillsat = setdiff(saturations, lower(satVars));
            nFill = numel(fillsat);
            assert(nFill == 0 || nFill == 1)
            if nFill == 1
                % Fill component is whichever saturation is assumed to fill
                % up the rest of the pores. This is done by setting that
                % increment equal to the negation of all others so that
                % sum(s) == 0 at end of update
                fillsat = fillsat{1};
                solvedFor = ~strcmpi(saturations, fillsat);
            else
                % All saturations are primary variables. Sum of saturations is
                % assumed to be enforced from the equation setup
                solvedFor = true(numel(saturations), 1);
            end
            ds = zeros(sum(state.nDof), numel(saturations));
            
            tmp = 0;
            ix = model.discretization.getDofIx(state, Inf);
            for phNo = 1:numel(saturations)
                if solvedFor(phNo)
                    v = model.getIncrement(dx, problem, saturations{phNo});
                    ds(ix, phNo) = v;
                    if nFill > 0
                        % Saturations added for active variables must be subtracted
                        % from the last phase
                        tmp = tmp - v;
                    end
                end
            end
            ds(ix, ~solvedFor) = tmp;
            % Use a more restrictive dsMaxAbs based on number of dofs
            dsAbsMax = model.parentModel.dsMaxAbs/min(model.discretization.basis.nDof, model.dsMaxAbsDiv);
            % We update all saturations simultanously, since this does not bias the
            % increment towards one phase in particular.
            state = model.updateStateFromIncrement(state, ds, problem, 'sdof', Inf, dsAbsMax);
            
        end
        
        % ----------------------------------------------------------------%
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            % Update state after convergence
            state.FacilityFluxProps = state.wellStateDG.FacilityFluxProps;
            % Let transport model do it's thing
            model.parentModel.outputFluxes = false;
            [state, report] = updateAfterConvergence@TransportModel(model, state0, state, dt, drivingForces);
            % Remove cell/face/well state
            state = rmfield(state, 'cellStateDG');
            state = rmfield(state, 'faceStateDG');
            state = rmfield(state, 'wellStateDG');
            % Limit solution
            if ~isempty(model.limiters)
                if model.storeUnlimited
                    % Store unlimited state if requested
                    state.ul = state;
                    if isfield(state.ul, 'ul')
                        state.ul = rmfield(state.ul, 'ul');
                    end
                end
                % Apply limiters
                for l = 1:numel(model.limiters)
                    limiter = model.limiters(l);
                    for v = 1:numel(limiter.variables)
                        state = limiter.function(state, limiter.variables{v}, limiter.tol, limiter.limits{v});
                    end
                end
            end
        end
        
        % ----------------------------------------------------------------%
        function dt = getMaximumTimestep(model, state, state0, dt0, drivingForces)
            % Define the maximum allowable time-step based on physics or
            % discretization choice
            state = model.getStateAD(state, false);
            state = state.wellStateDG;
            state.cells = (1:model.G.cells.num)';
            state.faces = (1:model.G.faces.num)';
            state.faces = state.faces(model.parentModel.operators.internalConn);
            dt = model.parentModel.getMaximumTimestep(state, state0, dt0, drivingForces);
        end
        
    end
    
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

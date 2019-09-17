classdef TransportModelDG < TransportModel
    
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function model = TransportModelDG(parent, varargin)
           
            model = model@TransportModel(parent);
            model.disc = [];
            [model, discArgs] = merge_options(model, varargin{:});
            % Construct discretization
            if isempty(model.disc)
                model.disc = DGDiscretization(model, discArgs{:});
            end
            model.parentModel.disc = model.disc;
            
            model.parentModel.operators = setupOperatorsDG(model.disc, model.parentModel.G, model.parentModel.rock);
            model.parentModel.outputFluxes = false;
            
        end
        
        %-----------------------------------------------------------------%
        function id = isDof(model, name)
            id = numel(name) > 3 && strcmp(name(end-2:end), 'dof');
        end
        
        
        
        %-----------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name, varargin)
            switch(lower(name))
                case {'swdof'}
                    index = model.satVarIndex('sw');
                    fn = 'sdof';
                case {'sodof'}
                    index = model.satVarIndex('so');
                    fn = 'sdof';
                case {'sgdof'}
                    index = model.satVarIndex('sg');
                    fn = 'sdof';
                case {'sdof'}
                    index = ':';
                    fn = 'sdof';
                case {'stdof'}
                    index = ':';
                    fn = 'sTdof';
                case {'pressuredof'}
                    index = ':';
                    fn = 'pressuredof';
                otherwise
                    [fn, index] = getVariableField@TransportModel(model, name, varargin{:});
            end
%             if ~isempty(fn)
%                 fn = [fn, 'dof'];
%             end
        end
        
        %-----------------------------------------------------------------%
        function index = satVarIndex(model, name)
            index = model.parentModel.satVarIndex(name);
        end
        
        %-----------------------------------------------------------------%
        function state = validateState(model, state)
            
            state.degree = repmat(model.disc.degree, model.G.cells.num, 1);
            wm = model.parentModel.FacilityModel.WellModels;
            for i = 1:numel(wm)
                state.degree(wm{i}.W.cells) = 0;
            end
            state = validateState@TransportModel(model, state);
            state = assignDofFromState(model.disc, state);
%             state = model.disc.getCellMean(state, state.sT);
            % TODO: Validate DG state props
        end
        
        %-----------------------------------------------------------------%
        function [state, names, origin] = getStateAD(model, state, init)
            if nargin < 3
                init = true;
            end
            parent = model.parentModel;
            % Get the AD state for this model
            [~, basenames, origin] = model.getPrimaryVariables(state);
            isParent = strcmp(origin, class(parent));
            
            basenames0 = basenames(isParent);
            basenames = cellfun(@(bn) [bn, 'dof'], basenames0, 'UniformOutput', false);
            origin = origin(isParent);
            basevars = cell(1, numel(basenames));
            for bNo = 1:numel(basenames)
                basevars{bNo} = model.getProp(state, basenames{bNo});
            end
            % Find saturations
            isS = false(size(basevars));
            nph = parent.getNumberOfPhases();
            phase_variable_index = zeros(nph, 1);
            for i = 1:numel(basevars)
                [f, ix] = model.getVariableField(basenames{i});
                if strcmp(f, 'sdof')
                    isS(i) = true;
                    phase_variable_index(ix) = i;
                end
            end
            % Figure out saturation logic
            isP = strcmp(basenames, 'pressuredof');
            vars = basevars;
            names = basenames;
            useTotalSaturation = strcmpi(model.formulation, 'totalSaturation') ...
                                    && sum(isS) == nph - 1;
            assert(useTotalSaturation, 'DG currently only supports total saturation formulation!');
            if useTotalSaturation
                % Replace pressure with total saturation
                replacement = 'sTdof';
                sTdof = model.getProp(state, replacement);
                % Replacing
                vars{isP} = sTdof;
                names{isP} = replacement;
                origin{isP} = class(model);
            else
                % Remove pressure and skip saturation closure
                vars = vars(~isP);
                names = names(~isP);
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
            % Let parent model handle state initialization
            state = model.initStateAD(state, basevars, basenames0, origin);
            
            for i = 1:numel(basenames)
                if any(strcmpi(basenames0{i}, {'sw', 'so', 'sg'}));
                   basenames0{i} = 's';
                   basenames{i}  = 'sdof';
                end
                v     = model.getProp(state, basenames0{i});
                state = model.setProp(state, basenames{i}, v);
                vm    = model.disc.getCellMean(state, value(v));
                state = model.setProp(state, basenames0{i}, vm); 
            end
            if useTotalSaturation
                % Set total saturation as well
                sTdof = vars{isP};
                state.sTdof = sTdof;
                % Evaluate at cell cubature points
                cellValue = model.disc.evaluateProp(state, sTdof, 'cell');
                state.cellStateDG.sT = cellValue;
                state.wellStateDG.sT = cellValue;
                % Evaluate at face cubature points
                faceValue = model.disc.evaluateProp(state, sTdof, 'face');
                state.faceStateDG.sT = faceValue;
                state.sT = model.disc.getCellMean(state, value(state.sTdof));
            end
        end
        
        function state = initStateAD(model, state, vars, names, origin)
              
            model.parentModel.G.cells.num = sum(state.nDof);
            state = initStateAD@TransportModel(model, state, vars, names, origin);
            state.sdof = state.s;
            state = model.evaluateBaseVariables(state);
            
        end
        
        function state = evaluateBaseVariables(model, state)
            
            [cellStateDG, faceStateDG, wellStateDG] = deal(state);
            
            names = {'pressure', 's'};
            for k = 1:numel(names)
                name = names{k};
                if isfield(state, name)
                    % Get dofs
                    dof = model.getProp(state, [name, 'dof']);
                    % Evaluate at cell cubature points
                    cellValue = model.disc.evaluateProp(state, dof, 'cell');
                    cellStateDG = model.setProp(cellStateDG, name, cellValue);
                    wellStateDG = model.setProp(wellStateDG, name, cellValue);
                    % Evaluate at face cubature points
                    faceValue = model.disc.evaluateProp(state, dof, 'face');
                    faceStateDG = model.setProp(faceStateDG, name, faceValue);
                end
            end
            
            cellStateDG.sT = getTotalSaturation(cellStateDG.s);
            wellStateDG.sT = getTotalSaturation(wellStateDG.s);
            faceStateDG.sT = getTotalSaturation(faceStateDG.s);
            
            
            [~, ~, cells] = model.disc.getCubature((1:model.G.cells.num)', 'volume');
            [~, ~, ~, faces] = model.disc.getCubature(find(model.parentModel.operators.internalConn), 'face');
            fcells = [model.G.faces.neighbors(faces,1); model.G.faces.neighbors(faces,2)];
            cellStateDG.type  = 'cell';
            cellStateDG.cells = cells;
            cellStateDG.fcells = fcells;
            cellStateDG.faces = faces;
            wellStateDG.type  = 'cell';
            cellStateDG.cells = cells;
            faceStateDG.type  = 'face';
            faceStateDG.cells = fcells;
            faceStateDG.faces = faces;
            
            state.cellStateDG = cellStateDG;
            state.wellStateDG = wellStateDG;
            state.faceStateDG = faceStateDG;
            
        end
        
        function state = assignBaseVariables(model, state)
            
            names = {'s'};
            for name = names
                if isfield(state, name{1}) && isfield(state, [name{1}, 'dof'])
                    dof = model.getProp(state, [name{1}, 'dof']);
                    v   = model.disc.getCellMean(state, dof);
                    state.(name{1}) = v;
                end
            end
             
            if strcmpi(model.formulation, 'totalSaturation')
                if isfield(state, 'sT') && isfield(state, 'sTdof')
                    dof = model.getProp(state, 'stdof');
                    v   = model.disc.getCellMean(state, dof);
                    state.sT = v;
                end
            end
             
        end
        
        function model = validateModel(model, varargin)
            model = validateModel@TransportModel(model, varargin{:});
                        
            model.parentModel.FluxDiscretization = FluxDiscretizationDG(model.parentModel);
            fp = model.parentModel.FlowPropertyFunctions;
            pvt = fp.getRegionPVT(model.parentModel);
            fp = fp.setStateFunction('PoreVolume', MultipliedPoreVolumeDG(model.parentModel, pvt));
            fp = fp.setStateFunction('GravityPermeabilityGradient', GravityPermeabilityGradientDG(model.parentModel));
            model.parentModel.FlowPropertyFunctions = fp;
            
        end
        
        %-----------------------------------------------------------------%
        function [eqs, names, types, state] = getModelEquations(tmodel, state0, state, dt, drivingForces)
            state0 = tmodel.evaluateBaseVariables(state0);
            model = tmodel.parentModel;
            [acc, flux, cellflux, names, types] = model.FluxDiscretization.componentConservationEquations(model, state, state0, dt);
            src = model.FacilityModel.getComponentSources(state);
            % Assemble equations and add in sources
            if strcmpi(tmodel.formulation, 'missingPhase')
                % Skip the last phase! Only mass-conservative for
                % incompressible problems
                acc = acc(1:end-1);
                flux = flux(1:end-1);
                names = names(1:end-1);
                types = types(1:end-1);
            end
            d        = tmodel.disc;
            d.nDof   = state.nDof;
            d.dofPos = state.dofPos;
            psi = d.basis.psi;
            grad_psi = d.basis.grad_psi;
            ix    = d.getDofIx(state, 1, src.cells);
            cells = rldecode((1:model.G.cells.num)', d.nDof, 1);
            d.sample = acc{1}(d.getDofIx(state, Inf));
            eqs = cell(1, numel(acc));
            for i = 1:numel(acc)
                eqs{i} = d.inner(acc{i}     , psi     , 'dV') ...
                       - d.inner(cellflux{i}, grad_psi, 'dV') ...
                       + d.inner(flux{i}    , psi     , 'dS');
                if ~isempty(src.cells)
                    eqs{i}(ix) = eqs{i}(ix) - src.value{i};
                end
                if ~model.useCNVConvergence
                    pv     = model.operators.pv(cells);
                    eqs{i} = eqs{i}.*(dt./pv);
                end    
            end
        end
        
        %-----------------------------------------------------------------%
        function [model, state] = prepareTimestep(model, state, state0, dt, drivingForces)
            [model, state] = prepareTimestep@TransportModel(model, state, state0, dt, drivingForces);
            state = assignDofFromState(model.disc, state, {'pressure'});
        end
        
        %-----------------------------------------------------------------%
        function [restVars, satVars, wellVars] = splitPrimaryVariables(model, vars)
            vars = cellfun(@(n) n(1:end-3), vars, 'UniformOutput', false);
            [restVars, satVars, wellVars] = model.parentModel.splitPrimaryVariables(vars);
            restVars = cellfun(@(n) [n, 'dof'], restVars, 'UniformOutput', false);
            satVars = cellfun(@(n) [n, 'dof'], satVars, 'UniformOutput', false);
        end
        
        %-----------------------------------------------------------------%
        function [state, report] = updateState(model, state, problem, dx, drivingForces)

            state = rmfield(state, 'cellStateDG');
            state = rmfield(state, 'faceStateDG');
            state = rmfield(state, 'wellStateDG');

            [restVars, satVars] = model.splitPrimaryVariables(problem.primaryVariables);
            
            % Update saturation dofs
            state = model.updateSaturations(state, dx, problem, satVars);
            % Update remaining non-saturation dofs
            state = model.updateDofs(state, dx, problem, restVars);
            % Update cell averages from dofs
            state0 = state;
            state  = model.assignBaseVariables(state);
            
            dx0 = model.getMeanIncrement(state, state0, problem);

            problem0 = problem;
            problem0.primaryVariables = cellfun(@(n) n(1:end-3), problem0.primaryVariables, 'UniformOutput', false);
            [state_tmp, report] = updateState@TransportModel(model, state0, problem0, dx0, drivingForces);
            dx_tmp = model.getMeanIncrement(state_tmp, state0, problem);
            
            frac = cellfun(@(x,y) x./y, dx_tmp, dx0, 'UniformOutput', false);
            cells = rldecode((1:model.G.cells.num)', state.nDof, 1);
            
%             vars = problem.primaryVariables;
%             for i = 1:numel(vars)
%                 f = frac{i};
%                 f(isnan(f)) = 1;
%                 f = f(cells);
%                 s     = model.getProp(state, vars{i});
%                 state = model.setProp(state, vars{i}, s.*f);
%             end
            
            
        end
        
        %-----------------------------------------------------------------%
        function dx = getMeanIncrement(model, state, state0, problem)
            
            vars = problem.primaryVariables;
            dx   = cell(1, numel(vars));
            for i = 1:numel(vars)
                vn = vars{i}(1:end-3);
                v  = model.getProp(state, vn);
                v0 = model.getProp(state0, vn);
                dx{i} = v - v0;
            end

        end
        
        % ----------------------------------------------------------------%
        function state = updateDofs(model, state, dx, problem, dofVars)
            
            for i = 1:numel(dofVars)
                state = updateStateFromIncrement(model, state, dx{i}, problem, dofVars{i}, inf, inf);
            end
            
        end
        
        % ----------------------------------------------------------------%
        function state = updateSaturations(model, state, dx, problem, satVars)

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
%             active = ~model.G.cells.ghost;
            ix = model.disc.getDofIx(state, Inf);
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
            % We update all saturations simultanously, since this does not bias the
            % increment towards one phase in particular.
            state   = model.updateStateFromIncrement(state, ds, problem, 'sdof', Inf, Inf);
            
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@TransportModel(model, state0, state, dt, drivingForces);
            state = rmfield(state, 'cellStateDG');
            state = rmfield(state, 'faceStateDG');
            state = rmfield(state, 'wellStateDG');
            
            propfn = model.parentModel.getStateFunctionGroupings();
            d = model.disc;
            d.nDof = state.nDof;
            d.dofPos = state.dofPos;
            ix = d.getDofIx(state, 1, Inf);
            psi    = model.disc.basis.psi(1);
            d.sample = state.sdof(:,1);
            for i = 1:numel(propfn)
                p = propfn{i};
                struct_name = p.getStateFunctionContainerName();
                names = p.getNamesOfStateFunctions();
                if isfield(state, struct_name)
                    for j = 1:numel(names)
                        name = names{j};
                        if ~isempty(state.(struct_name).(name))
                            v = state.(struct_name).(name);
                            nph = numel(v);
                            for ph = 1:nph
                                v{ph} = d.inner(v{ph}, psi, 'dV');
                                v{ph} = v{ph}(ix);
                            end
                            state.(struct_name).(name) = v;
                        end
                    end
                end
            end
            
        end
        
    end
    
end

function sT = getTotalSaturation(s)
    if iscell(s)
        sT  = 0;
        nph = numel(s);
        for i = 1:nph
            sT = sT + s{i};
        end
    else
        sT = sum(s,2);
    end
end
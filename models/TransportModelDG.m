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
            
        end
        
        %-----------------------------------------------------------------%
        function id = isDof(model, name)
            id = numel(name) > 3 && strcmp(name(end-2:end), 'dof');
        end
        
        %-----------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name, varargin)
            [fn, index] = getVariableField@TransportModel(model, name, varargin{:});
            if ~isempty(fn)
                fn = [fn, 'dof'];
            end
        end
        
        %-----------------------------------------------------------------%
        function state = validateState(model, state)
            
            state.degree = repmat(model.disc.degree, model.G.cells.num, 1);
            wm = model.parentModel.FacilityModel.WellModels;
            for i = 1:numel(wm)
                state.degree(wm{i}.W.cells) = 0;
            end
            state    = assignDofFromState(model.disc, state);
            state    = validateState@TransportModel(model, state);
            state.sT = model.disc.getCellMean(state, state.sTdof);
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
            basenames = basenames(isParent);
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
            isP = strcmp(basenames, 'pressure');
            vars = basevars;
            names = basenames;
            useTotalSaturation = strcmpi(model.formulation, 'totalSaturation') ...
                                    && sum(isS) == nph - 1;
            assert(useTotalSaturation, 'DG currently only supports total saturation formulation!');
            if useTotalSaturation
                % Replace pressure with total saturation
                replacement = 'sT';
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
            state = model.initStateAD(state, basevars, basenames, origin);
            state.s = model.disc.getCellMean(state, value(state.sdof));
            state.pressure = model.disc.getCellMean(state, value(state.pressuredof));
            if useTotalSaturation
                % Set total saturation as well
                sTdof = vars{isP};
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
                if isfield(state, [name, 'dof'])
                    % Get dofs
                    dof = model.getProp(state, name);
                    % Evaluate at cell cubature points
                    cellValue = model.disc.evaluateProp(state, dof, 'cell');
                    cellStateDG = model.parentModel.setProp(cellStateDG, name, cellValue);
                    wellStateDG = model.parentModel.setProp(wellStateDG, name, cellValue);
                    % Evaluate at face cubature points
                    faceValue = model.disc.evaluateProp(state, dof, 'face');
                    faceStateDG = model.parentModel.setProp(faceStateDG, name, faceValue);
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
        function [acc, names, types, state] = getModelEquations(tmodel, state0, state, dt, drivingForces)
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
                if ~isempty(src.cells)
                    acc{i}(ix) = acc{i}(ix) - src.value{i};
                end
                eqs{i} = d.inner(acc{i}     , psi     , 'dV') ...
                       - d.inner(cellflux{i}, grad_psi, 'dV') ...
                       + d.inner(flux{i}    , psi     , 'dS');
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
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            state_dof   = state;
            state_dof.s = state.sdof;
            model.parentModel.G.cells.num = sum(state.nDof);
            [state_dof, report] = updateState@TransportModel(model, state_dof, problem, dx, drivingForces);
            state.sdof = state_dof.s;
            state.s = model.disc.getCellMean(state, state.sdof);
            state = rmfield(state, 'cellStateDG');
            state = rmfield(state, 'faceStateDG');
            state = rmfield(state, 'wellStateDG');
        end
        
        % ----------------------------------------------------------------%
        function state = updateSaturations(model, state, dx, problem, satDofVars)

            if nargin < 5
                % Get the saturation names directly from the problem
                [~, satDofVars] = ...
                    splitPrimaryVariables(model, problem.primaryVariables);
            end
            if isempty(satDofVars)
                % No saturations passed, nothing to do here.
                return
            end
            % Solution variables should be saturations directly, find the missing
            % link
            saturations = lower(model.getDGDofVarNames);
            
            fillsat = setdiff(saturations, lower(satDofVars));
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
            active = ~model.G.cells.ghost;
            ix = model.disc.getDofIx(state, Inf, active);
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
            state   = model.updateStateFromIncrement(state, ds, problem, 'sdof', Inf, model.dsMaxAbs);
            state.s = model.disc.getCellSaturation(state);

%             ix = any(abs(ds)>model.dsMaxAbs,2);
%             alph = model.dsMaxAbs./max(abs(ds(ix,:)), [], 2);
%             
%             ds(ix,:) = alph.*ds(ix,:);
%             state   = model.updateStateFromIncrement(state, ds, problem, 'sdof', Inf, Inf);
%             state.s = model.disc.getCellSaturation(state);            
            
            if nFill == 1
                
                if 1
                bad = any((state.s > 1 + model.disc.meanTolerance) ...
                        | (state.s < 0 - model.disc.meanTolerance), 2);
                
                else
                [smin, smax] = model.disc.getMinMaxSaturation(state);
                over  = smax > 1 + model.disc.meanTolerance;
                under = smin < 0 - model.disc.meanTolerance;
                bad = over | under;
                end
                    
                    
                if any(bad)
                    state.s(bad, :) = min(state.s(bad, :), 1);
                    state.s(bad, :) = max(state.s(bad, :), 0);
                    state.s(bad, :) = bsxfun(@rdivide, state.s(bad, :), ...
                                                  sum(state.s(bad, :), 2));
                    state = dgLimiter(model.disc, state, bad, 's', 'kill');
                end
            else
                bad = any(state.s < 0 - model.disc.meanTolerance, 2);
                 if any(bad)
                    state.s(bad, :) = max(state.s(bad, :), 0);
                    state = dgLimiter(model.disc, state, bad, 's', 'kill');
                 end
            end

            if model.disc.limitAfterNewtonStep
                % Limit solution
                state = model.disc.limiter(model, state, [], true);
            end
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@TransportModel(model, state0, state, dt, drivingForces);
            state = rmfield(state, 'cellStateDG');
            state = rmfield(state, 'faceStateDG');
            state = rmfield(state, 'wellStateDG');
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
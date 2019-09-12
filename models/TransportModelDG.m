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
            state    = assignDofFromState(model.disc, state);    
            state    = validateState@TransportModel(model, state);
            state = rmfield(state, 'sTdof');
            state.sT = sum(state.s,2);
            state    = assignDofFromState(model.disc, state);
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
                basevars{bNo} = model.getProp(state, basenames{bNo}, false);
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
            if useTotalSaturation
                % Replace pressure with total saturation
                replacement = 'sT';
                sTdof = state.sTdof;
%                 sT = model.getProp(state, [replacement, 'dof'], false);
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
            state.sdof = state.s;
            state.s = model.disc.getCellMean(state, value(state.sdof));
            if useTotalSaturation
                % Set total saturation as well
                sT = vars{isP};
                state = model.setProp(state, replacement, sT);
            end
        end
        
        function model = validateModel(model, varargin)
            defaultedDiscretization = isempty(model.parentModel.FluxDiscretization);
            model = validateModel@TransportModel(model, varargin{:});
            if defaultedDiscretization
                pmodel = model.parentModel;
                fd = pmodel.FluxDiscretization;
                % Set flow state builder
                fd = fd.setFlowStateBuilder(FlowStateBuilderDG());
                % Replace some existing properties with DG variants
                fd = fd.setStateFunction('Pressure', PrimaryVariableDG(fd.Pressure));
                fd = fd.setStateFunction('GravityPotentialDifference', GravityPotentialDifferenceDG(pmodel));
                fd = fd.setStateFunction('TotalFlux', FixedTotalFluxDG(pmodel));
                
                
                fp = pmodel.FlowPropertyFunctions;
                % Replace some existing properties with DG variants
                fp = fp.setStateFunction('PhaseSaturations', PrimaryVariableDG(fp.PhaseSaturations));
                fp = fp.setStateFunction('Pressure', PrimaryVariableDG(fp.Pressure));
                
                % Replace object
                model.parentModel.FluxDiscretization    = fd;
                model.parentModel.FlowPropertyFunctions = fp;
            end
        end
        
        %-----------------------------------------------------------------%
        function [eqs, names, types, state] = getModelEquations(tmodel, state0, state, dt, drivingForces)
            model = tmodel.parentModel;
            [eqs, flux, names, types] = model.FluxDiscretization.componentConservationEquations(model, state, state0, dt);
            src = model.FacilityModel.getComponentSources(state);
            % Assemble equations and add in sources
            if strcmpi(tmodel.formulation, 'missingPhase')
                % Skip the last phase! Only mass-conservative for
                % incompressible problems
                eqs = eqs(1:end-1);
                flux = flux(1:end-1);
                names = names(1:end-1);
                types = types(1:end-1);
            end
            for i = 1:numel(eqs)
                if ~isempty(src.cells)
                    eqs{i}(src.cells) = eqs{i}(src.cells) - src.value{i};
                end
                eqs{i} = model.operators.AccDiv(eqs{i}, flux{i});
                 if ~model.useCNVConvergence
                    pv     = model.operators.pv;
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
            [state_dof, report] = updateState@TransportModel(model, state_dof, problem, dx, drivingForces);
            state.sdof = state_dof.s;
            state.s = model.disc.getCellMean(state, state.sdof);
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
        
    end
    
end

% --------------------------------------------------------------------%
function [fn, index] = getVariableField(model, name, varargin)

    if strcmpi(name(end-2:end), 'dof')
        nm = name(1:end-3);
    end
    [fn, index] = model.getVariableField(nm, varargin{:});
    fn = [fn, 'dof'];

end
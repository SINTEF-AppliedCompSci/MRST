classdef TransportModel < WrapperModel
    properties
        formulation = 'totalSaturation'
    end
    
    methods
        function model = TransportModel(parent, varargin)
            model = model@WrapperModel(parent);
            model = merge_options(model, varargin{:});
            model.AutoDiffBackend = parent.AutoDiffBackend;
        end
        
        function [state, names, origin] = getStateAD(model, state, init)
            if nargin < 3
                init = true;
            end
            parent = model.parentModel;
            % Get the AD state for this model
            [basevars, basenames, baseorigin] = model.getPrimaryVariables(state);
            isParent = strcmp(baseorigin, class(parent));
            basevars = basevars(isParent);
            basenames = basenames(isParent);
            baseorigin = baseorigin(isParent);
            % Find saturations
            isS = false(size(basevars));
            nph = parent.getNumberOfPhases();
            phase_variable_index = zeros(nph, 1);
            for i = 1:numel(basevars)
                bn = basenames{i};
                if strcmp(bn, 'x')
                    % Hack for blackoil-models
                    ss = true;
                else
                    [f, ix] = model.getVariableField(basenames{i});
                    ss = strcmp(f, 's');
                end
                if ss
                    isS(i) = true;
                    phase_variable_index(ix) = i;
                end
            end
            % Figure out saturation logic
            isP = strcmp(basenames, 'pressure');
            vars = basevars;
            names = basenames;
            origin = baseorigin;
            useTotalSaturation = strcmpi(model.formulation, 'totalSaturation') ...
                                    && sum(isS) == nph - 1;
            if useTotalSaturation
                % Replace pressure with total saturation
                replacement = 'sT';
                sT = model.getProp(state, replacement);
                % Replacing
                vars{isP} = sT;
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
            state = model.initStateAD(state, basevars, basenames, baseorigin);
            if useTotalSaturation
                % Set total saturation as well
                sT = vars{isP};
                state = model.setProp(state, replacement, sT);
            end
        end
        

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

        function state = validateState(model, state)
            state = validateState@WrapperModel(model, state);
            fn = model.getVariableField('sT');
            if isfield(state, 's') && ~isfield(state, fn)
                state = model.setProp(state, 'sT', sum(state.s, 2));
            end
        end
        
        function model = validateModel(model, varargin)
            defaultedDiscretization = isempty(model.parentModel.FluxDiscretization);
            model = validateModel@WrapperModel(model, varargin{:});
            if defaultedDiscretization
                pmodel = model.parentModel;
                fd = pmodel.FluxDiscretization;
                fp = pmodel.FlowPropertyFunctions;
                % Replace existing properties with total flux variants
                fd = fd.setStateFunction('PhaseFlux', PhaseFluxFixedTotalVelocity(pmodel));
                fd = fd.setStateFunction('PhaseUpwindFlag', PhasePotentialUpwindFlag(pmodel));
                fd = fd.setStateFunction('ComponentPhaseFlux', ComponentPhaseFluxFractionalFlow(pmodel));
                % Set extra props
                fd = fd.setStateFunction('PhaseInterfacePressureDifferences', PhaseInterfacePressureDifferences(pmodel));
                fd = fd.setStateFunction('TotalFlux', FixedTotalFlux(pmodel));
                fd = fd.setStateFunction('FaceTotalMobility', FaceTotalMobility(pmodel));
                % Set flow properties
                if strcmpi(model.formulation, 'totalSaturation')
                    fp = fp.setStateFunction('TotalSaturation', TotalSaturation(pmodel));
                    fp = fp.setStateFunction('ComponentMobility', ComponentMobilityTotalSaturation(pmodel));
                    fp = fp.setStateFunction('ComponentPhaseMass', ComponentPhaseMassTotalSaturation(pmodel));
                end
                % Replace object
                model.parentModel.FluxDiscretization = fd;
                model.parentModel.FlowPropertyFunctions = fp;
            end
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = getEquations@PhysicalModel(model, state0, state, dt, drivingForces, varargin{:});
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            isS = strcmpi(problem.primaryVariables, 'sT');
            if any(isS)
                state = model.updateStateFromIncrement(state, dx, problem, 'sT', inf, inf);
                state = model.capProperty(state, 'sT', 1e-8);
                dx = dx(~isS);
            end
            problem.primaryVariables = problem.primaryVariables(~isS);
            
            [state, report] = model.parentModel.updateState(state, problem, dx, drivingForces);
        end
        
        function  [convergence, values, names] = checkConvergence(model, problem)
            [convergence, values, names] = model.parentModel.checkConvergence(problem);
        end
        
        function [fn, index] = getVariableField(model, name, varargin)
            if strcmpi(name, 'sT')
                fn = 'sT';
                index = ':';
            else
                [fn, index] = model.parentModel.getVariableField(name, varargin{:});
            end
        end
        
        function [model, state] = prepareTimestep(model, state, state0, dt, drivingForces)
            % Prepare state and model (temporarily) before solving a time-step
            [model, state] = prepareTimestep@WrapperModel(model, state, state0, dt, drivingForces);
        end

        function [state, report] = updateAfterConvergence(model, varargin)
            [state, report] = updateAfterConvergence@WrapperModel(model, varargin{:});
        end
    end
end
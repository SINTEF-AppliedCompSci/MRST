classdef GenericBlackOil < ThreePhaseBlackOilModel & ExtendedReservoirModel
    properties
        
    end
    
    methods
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            opt = struct('Verbose',     mrstVerbose,...
                        'reverseMode', false,...
                        'resOnly',     false,...
                        'iteration',   -1, ...
                        'drivingForces0', []);
            opt = merge_options(opt, varargin{:});


            % Define primary variables
            if opt.reverseMode
                [state0, primaryVars] = model.getReverseStateAD(state0);
                % The model must be validated with drivingForces so that the
                % FacilityModel gets updated.
                model = model.validateModel(drivingForces);
                state = model.getStateAD(state, false);
            else
                [state, primaryVars] = model.getStateAD(state, ~opt.resOnly);
            end
            [eqs, names, types, state] = model.getModelEquations(state0, state, dt, drivingForces);

            state = value(state);
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
        end
        
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces)
            [eqs, divTerms, names, types] = model.FluxDiscretization.componentConservationEquations(model, state, state0, dt);

            
            src = model.FacilityModel.getComponentSources(state);
            for i = 1:numel(divTerms)
                eqs{i} = eqs{i} + divTerms{i};
                eqs{i}(src.cells) = eqs{i}(src.cells) + src.value{i};
            end
            [weqs, wnames, wtypes, state] = model.FacilityModel.getModelEquations(state0, state, dt, drivingForces);
            
            rhoS = model.getSurfaceDensities();
            for i = 1:numel(divTerms)
                eqs{i} = eqs{i}./rhoS(i);
            end
            eqs = [eqs, weqs];
            names = [names, wnames];
            types = [types, wtypes];
        end
        
        function names = getComponentNames(model)
            names = cellfun(@(x) x.name, model.Components, 'UniformOutput', false);
        end
        
        function n = getNumberOfComponents(model)
            n = numel(model.Components);
        end
        
        function n = getNumberOfPhases(model)
            n = model.water + model.oil + model.gas;
        end
    end
end
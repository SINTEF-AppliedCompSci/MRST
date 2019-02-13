classdef GenericBlackOil < ThreePhaseBlackOilModel & ExtendedReservoirModel
    properties
        
    end
    
    methods
    function model = GenericBlackOil(G, rock, fluid, varargin)
        model = model@ThreePhaseBlackOilModel(G, rock, fluid, varargin{:});
        model.OutputProperties = {'ComponentTotalMass'};
        
        nph = model.getNumberOfPhases();
        model.Components = cell(1, nph);
        names = model.getPhaseNames();
        for ph = 1:nph
            switch names(ph)
                case 'W'
                    c = ImmiscibleComponent('water', ph);
                case 'O'
                    c = OilComponent('oil', ph);
                case 'G'
                    c = GasComponent('gas', ph);
                otherwise
                    error('Unknown phase');
            end
            model.Components{ph} = c;
        end
    end
        
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
            [eqs, flux, names, types] = model.FluxDiscretization.componentConservationEquations(model, state, state0, dt);
            src = model.FacilityModel.getComponentSources(state);
            % Assemble equations and add in sources
            for i = 1:numel(eqs)
                if ~isempty(src.cells)
                    eqs{i}(src.cells) = eqs{i}(src.cells) - src.value{i};
                end
                eqs{i} = model.operators.AccDiv(eqs{i}, flux{i});
            end
            % Get facility equations
            [weqs, wnames, wtypes, state] = model.FacilityModel.getModelEquations(state0, state, dt, drivingForces);
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
        
        function [state, report] = updateState(model, state, problem, dx, forces) 
            [state, report] = updateState@ThreePhaseBlackOilModel(model, state, problem, dx, forces);
            if ~isempty(model.FacilityModel)
                state = model.FacilityModel.applyWellLimits(state);
            end
        end
        
        function model = validateModel(model, varargin)
            % Validate model.
            %
            % SEE ALSO:
            %   :meth:`ad_core.models.PhysicalModel.validateModel`
            model.FacilityModel = ExtendedFacilityModel(model);
            model = validateModel@ThreePhaseBlackOilModel(model, varargin{:});
        end
    end
end
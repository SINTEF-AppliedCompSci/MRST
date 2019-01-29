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

            
            src = model.FacilityModel.getProp(state, 'ComponentTotalFlux');
            for i = 1:numel(divTerms)
%                 src = model.FacilityModel.getComponentSourceTerm(state, names{i});
                eqs{i} = eqs{i} + divTerms{i} + src;
            end
            [weqs, names, types, state] = model.FacilityModel.getModelEquations(state0, state, dt, drivingForces);
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
        
        % Various hacked versions of functions to bootstrap coming changes
        function [eqs, names, types, wellSol, src] = insertWellEquations(model, eqs, names, ...
                                                         types, wellSol0, wellSol, ...
                                                         wellVars, wellMap, ...
                                                         p, mob, rho, ...
                                                         dissolved, components, ...
                                                         dt, opt)
            if model.FacilityModel.getNumberOfActiveWells(wellSol) == 0
                src = [];
                return
            end
            fm = model.FacilityModel;
            nPh = nnz(model.getActivePhases);
            [src, wellsys, wellSol] = ...
                fm.getWellContributions(wellSol0, wellSol, wellVars, ...
                                        wellMap, p, mob, rho, dissolved, components, ...
                                        dt, opt.iteration);

            rhoS = model.getSurfaceDensities();
            wc = src.sourceCells;
            [~, longNames] = getPhaseNames(model);
            % Treat phase pseudocomponent source terms from wells
            for i = 1:nPh
                sub = strcmpi(names, longNames{i});
                if any(sub)
                    assert(strcmpi(types{sub}, 'cell'), 'Unable to add source terms to equation that is not per cell.');
                    eqs{sub}(wc) = eqs{sub}(wc) - src.phaseMass{i};
                end
            end
            eqs = horzcat(eqs, wellsys.wellEquations, {wellsys.controlEquation});
            names = horzcat(names, wellsys.names, 'closureWells');
            types = horzcat(types, wellsys.types, 'well');
        end
    end
end
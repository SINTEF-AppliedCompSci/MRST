classdef BlackOilFluidModel < ThreePhaseBlackOilModel

    properties
        primaryVars
    end

    methods
        function model = BlackOilFluidModel(G, rock, fluid, varargin)
            model = model@ThreePhaseBlackOilModel(G, rock, fluid);
            fluidModelType = 'blackoil';
            model.disgas = true;
            model.vapoil = false;
            primaryVars = []; % updated by MechFluidModel, used by the updateState function
            model = merge_options(model, varargin{:});
        end

        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, ...
                                                        varargin)
            error(['The BlackOilModel is not meant for being called directly. ' ...
                   'It should be used coupled with a mechanical model, derived ' ...
                   'from MechFluidModel.'])
        end

        function fds = getAllVarsNames(model)
            fds = {'wellSol', 'pressure', 's', 'rs', 'rv', 'sW', 'sG', 'sO', ...
                   'water', 'oil', 'gas'};
        end
        
        function [state, report] = updateState(model, state, problem, dx, ...
                                               drivingForces)
            vars = problem.primaryVariables;
            removed = true(size(vars));
            [lia, loc] = ismember(primaryVars, vars);
            assert(all(lia), 'The primary variables are not set correctly.');
            removed(loc) = false;
            problem.primaryVariables = vars(removed);
            dx   = dx(removed);
            
            [state, report] = updateState@ThreePhaseBlackOilModel(model, state, ...
                                                              problem, dx, drivingForces);
            
        end


    end

end

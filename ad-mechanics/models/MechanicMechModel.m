classdef MechanicMechModel < MechanicModel
% Mechanical model to be used with fully coupled solver
    properties
        primaryVars;
    end

    methods
        function model = MechanicMechModel(G, rock, mech_problem, varargin)
            model = model@MechanicalModel(G, rock, mech_problem, varargin{:});
            primaryVars = [];
            model = merge_options(model, varargin{:});

        end

        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, varargin)
            error(['The MechanicslMechModel is not meant for being called directly. ' ...
                   'It should be used coupled with a mechanical model, derived ' ...
                   'from MechFluidModel.'])
        end

        function [primaryVars, fds] = getAllVarsNames(model)
            primaryVars = {'xd'};
            fds = {'xd', 'uu', 'u', 'stress', 'strain', 'vdiv'};
        end


        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            vars = problem.primaryVariables;
            removed = true(size(vars));
            [lia, loc] = ismember(primaryVars, vars);
            assert(all(lia), 'The primary variables are not set correctly.');
            removed(loc) = false;
            problem.primaryVariables = vars(removed);
            dx   = dx(removed);

            [state, report] = updateState@MechanicalModel(model, state, problem, ...
                                                          dx, drivingForces);
            state = addDerivedQuantities(model, state);
        end

    end
end

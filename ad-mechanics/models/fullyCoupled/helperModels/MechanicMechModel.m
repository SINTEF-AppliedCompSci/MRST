classdef MechanicMechModel < MechanicModel
% Mechanical model to be used with fully coupled solver
    properties
        primaryVarNames;
    end

    methods
        function model = MechanicMechModel(G, rock, mech_problem, varargin)
            model = model@MechanicModel(G, rock, mech_problem, varargin{:});
            model.primaryVarNames = {'xd'};
            model = merge_options(model, varargin{:});

        end

        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, varargin)
            error(['The MechanicslMechModel is not meant for being called directly. ' ...
                   'It should be used coupled with a mechanical model, derived ' ...
                   'from MechFluidModel.'])
        end

        function  fds = getAllVarsNames(model)
            fds = {'xd', 'uu', 'u', 'stress', 'strain', 'vdiv'};
        end

        function varnames = getAllPrimaryVariables(model)
            varnames = model.primaryVarNames;
        end

        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            vars = problem.primaryVariables;
            ind = false(size(vars));
            mechvars = model.getAllPrimaryVariables();
            [lia, loc] = ismember(mechvars, vars);
            assert(all(lia), 'The primary variables are not set correctly.');

            ind(loc) = true;
            problem.primaryVariables = vars(ind);
            dx   = dx(ind);

            [state, report] = updateState@MechanicModel(model, state, problem, ...
                                                        dx, drivingForces);
            
            state = addDerivedQuantities(model, state);
        end

    end
end

function state = initDisplacement(model, state, uu, varargin)

    opt = struct('mech_equil', true, ...
                 'pressure'  , []   );
    opt = merge_options(opt, varargin{:});

    if (opt.mech_equil)

        if isprop(model, 'mechModel')
            mechmodel = model.mechModel; 
        else
            mechmodel = MechanicBiotModel(model.G, model.rock, model.mech, ...
                                          'InputModel', model);
        end

        state = mechmodel.setProp(state, 'xd', mechmodel.operators.mech.rhs); % Dummy values, just used
                                                                              % to get the correct dimension.

        if ~isempty(opt.pressure)
            drivingForces.fluidp = opt.pressure;
        else
            drivingForces.fluidp = 0;
        end

        solver = NonLinearSolver(); % The problem is linear but we use the general
                                    % framework. The cost should be
                                    % minimal.
        [state, failure, report] = ...
            solveMinistep(solver, mechmodel, state, state, 0, drivingForces);

    else
        error('not checked')
        u = reshape(uu', [], 1);
        state.xd = zeros(size(model.operators.mech.A, 2), 1);
        state.xd = u(~model.operators.mech.isdirdofs);
    end

    state = addDerivedQuantities(mechmodel, state);

end

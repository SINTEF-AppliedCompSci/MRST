function state = computeInitDisp(model, state, uu, varargin)

    opt = struct('mech_equil', true, ...
                 'pressure'  , []   );
    opt = merge_options(opt, varargin{:});

    mechModel = model.mechModel;
    
    if (opt.mech_equil)

        mechmodel = MechanicModel(mechModel.G, mechModel.rock, mechModel.mech);

        state = mechmodel.setProp(state, 'xd', mechmodel.operators.rhs); % Dummy values, just used
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
        state.xd = zeros(size(model.operators.A, 2), 1);
        state.xd = u(~model.operators.isdirdofs);
    end

    state = addDerivedQuantities(mechmodel, state);

end

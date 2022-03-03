function nonLinearSolver = SetupNonLinearSolver()
%     params = model.solver.params;
    params = SolverParams();
    % Time Step Selector
    timestepper = StateChangeTimeStepSelector(...
                  'firstRampupStepRelative', params.firstRampupStepRelative,...
                  ...'firstRampupStep',         params.firstRampupStep,...
                  ...'targetIterationCount',    params.targetIterationCount,...
                  ...'iterationOffset',         params.iterationOffset,...
                  'targetProps',             params.targetProps,...
                  ...'targetChangeRel',         params.targetChangeRel,...
                  'targetChangeAbs',         params.targetChangeAbs);
                                           
    % linear solver
    linsolve = [];
    if( params.useCPR )
        if ~isempty(mrstPath('agmg'))
            mrstModule add agmg
            pSolver = AGMGSolverAD();
        else
            pSolver = BackslashSolverAD();
        end
        linsolve = CPRSolverAD('ellipticSolver', pSolver);
    elseif( params.useGMRES )
        linsolve = GMRES_ILUSolverAD();
    end
    
    % nonlinear solver
    nonLinearSolver = NonLinearSolver('LinearSolver', linsolve,...
                                      'timeStepSelector', timestepper,...
                                      'maxIterations',    params.maxIterations,...
                                      'minIterations',    params.minIterations,...
                                      'maxTimestepCuts',  params.maxTimestepCuts);
end

function params = SolverParams()
    params = struct( ...
        ... %% Linear Solver settings
        ... % Solve a problem with a pressure component using constrained a pressure residual method
        'useCPR',                   false,      ...
        ... % Preconditioned GMRES solver
        'useGMRES',                 false,      ...
        ... %% Nonlinear Solver settings
        ... % The maximum number of iterations during a ministep of nonlinear solver
        'maxIterations',            100,        ... % default 25
        ... % The minimum number of iterations during a ministep o fnonlinear solver
        'minIterations',            2,          ... % default 1 
        ... % The maximum number of times the timestep can be halved before it
        ... % is counted as a failed run
        'maxTimestepCuts',          100,        ... % default 6
        ... %% Time Step Selector settings
        ... % The target properties to time control. Cell array of N strings,
        ... % each of which are suitable for the simulation model's getProp
        ... % function. Typically, this means that the property name has been
        ... % implemented in the Model's "getVariableField" member function.
        ... % See PhysicalModel for more information.
        'targetProps',              {'s'},      ... % default {}
        ... % Target change (relative units). This is a double array of the
        ... % same size as targetProps, where the relative change in properties
        ... % are set to these targets. The i-th property will use the i-th
        ... % entry of targetChangeRel. Set values to inf to not apply relative
        ... % changes.
        'targetChangeRel',          0.01,       ... % default 1
        ... % Target change (absolute units). Double array of the same length
        ... % as targetProps. The absolute units are useful when dealing with
        ... % non-scaled values.
        'targetChangeAbs',          0.01,       ... % default inf
        ... % Desired number of nonlinear iterations per timestep.
        'targetIterationCount',     5,          ... % default 5 
        ... % Offset to make iteration a bit smoother as a response function.
        'iterationOffset',          5,          ... % default 5
        ... % The first ministep attempted after controls have changed. Could
        ... % be set to a low value to get the timestep controller started with
        ... % some estimate of problem stiffness.
        'firstRampupStep',          [],         ... % default inf
        ... % Same as firstRampupStep, but interpreted in a relative fashion
        ... % (i.e. if timestep is 5*days, the ramp up step will then be
        ... % 5*days*firstRampupStepRelative.
        'firstRampupStepRelative',  0.1         ... % default 1
    );
end
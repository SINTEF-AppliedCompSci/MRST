mrstModule add ad-core ad-blackoil spe10 blackoil-sequential mrst-gui linearsolvers

% Select layer 1
layers = 1;
mrstModule add ad-core ad-blackoil blackoil-sequential spe10

% The base case for the model is 2000 days. This can be reduced to make the
% problem faster to run.
T = 2000*day;
[state, model, schedule] = setupSPE10_AD('layers', layers, 'dt', 30*day, 'T',  T);

% Set up pressure and transport linear solvers
% AMGCL in AMG mode with default parameters
psolver = AMGCLSolverAD('coarsening', 'smoothed_aggregation', 'maxIterations', 50, 'tolerance', 1e-4, 'relaxation', 'ilu0');
psolver.keepNumber = model.G.cells.num;
% psolver.applyRightDiagonalScaling = true;
% AMGCL without AMG as a Krylov solver with ILU(0) preconditioner
tsolver = AMGCLSolverAD('preconditioner', 'relaxation', 'relaxation', 'ilu0', 'tolerance', 1e-4);

% Set up the sequential model
seqModel = getSequentialModelFromFI(model, ...
    'pressureLinearSolver', psolver,....
    'transportLinearSolver', tsolver);

% We set up a timestep selector that aims for timesteps where the
% maximum saturation change is equal to a fixed value.
stepSel = StateChangeTimeStepSelector(...
    'targetProps', {'s'},...
    'targetChangeAbs', 0.25);
% Run problem
solver = NonLinearSolver('timeStepSelector', stepSel);
[wsSeq, statesSeq, repSeq] = simulateScheduleAD(state, seqModel, schedule, 'NonLinearSolver', solver);

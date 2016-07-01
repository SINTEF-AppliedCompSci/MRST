%% Solve layers of SPE10 with sequential and fully implicit solver
% This is a short example demonstrating how to solve any number of layers
% of the SPE10, model 2 problem using both a sequential and a fully
% implicit solver.
mrstModule add ad-core ad-blackoil spe10 blackoil-sequential mrst-gui

% Set up pressure and transport linear solvers
if ~isempty(mrstPath('agmg'))
    mrstModule add agmg
    psolver = AGMGSolverAD();
else
    psolver = BackslashSolverAD();
end
tsolver = GMRES_ILUSolverAD();
% Select layer 1
layers = 1;
mrstModule add ad-core ad-blackoil blackoil-sequential spe10

% The base case for the model is 2000 days. This can be reduced to make the
% problem faster to run.
T = 2000*day;
[state, model, schedule] = setupSPE10_AD('layers', layers, 'dt', 30*day, ...
                                                           'T',  T);
% Set up the sequential model
seqModel = getSequentialModelFromFI(model, 'pressureLinearSolver', psolver,....
                                           'transportLinearSolver', tsolver);
% We set up a timestep selector that aims for timesteps where the
% maximum saturation change is equal to a fixed value.
stepSel = StateChangeTimeStepSelector('targetProps', {'s'},...
                                      'targetChangeAbs', 0.25);
% Run problem
solver = NonLinearSolver('timeStepSelector', stepSel);
[wsSeq, statesSeq, repSeq] = simulateScheduleAD(state, seqModel, schedule, 'NonLinearSolver', solver);
%% Simulate fully implicit
% Solve the fully implicit version of the problem, with a CPR
% preconditioner that uses the same linear solver for the pressure as the
% sequential solver.
solver.timeStepSelector.reset();
solver.LinearSolver = CPRSolverAD('ellipticSolver', psolver);

[wsFIMP, statesFIMP, repFIMP] = simulateScheduleAD(state, model, schedule, 'NonLinearSolver', solver);
%% Plot simulation time taken
% The sequential solver can be much faster than the fully implicit solver
% for certain problems.
figure(1); clf; hold on
plot(cumsum(repSeq.SimulationTime)/60, 'b-*')
plot(cumsum(repFIMP.SimulationTime)/60, 'r-o')
ylabel('Simulation time [minutes]')
xlabel('Control step #')
legend('Sequential implicit', 'Fully implicit')
%% Plot the results in interactive viewers
G = model.G;
W = schedule.control(1).W;

% Plot the well curves
plotWellSols({wsSeq, wsFIMP}, cumsum(schedule.step.val), ...
            'datasetnames', {'Sequential', 'FIMP'}, 'field', 'qOs')

% Plot reservoir quantities
figure;
plotToolbar(G, statesSeq);
axis equal tight
view(90, 90);
plotWell(G, W);
title('Sequential implicit')

figure;
plotToolbar(G, statesFIMP);
axis equal tight
view(90, 90);
plotWell(G, W);
title('Fully-implicit')
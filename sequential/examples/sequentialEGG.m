%% Simulate the Egg model
% The Egg model is a nearly incompressible waterflood problem. This example
% sets up a sequential solver and compares the results with the fully
% implicit discretization. As there are small changes in mobility, this
% problem is very well suited to sequential splitting.
%
% For details on the EggModel, see Jansen, J. D., et al. "The egg model–a
% geological ensemble for reservoir simulation." Geoscience Data Journal
% 1.2 (2014): 192-195.

mrstModule add ad-core ad-blackoil spe10 sequential mrst-gui example-suite

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
mrstModule add ad-core ad-blackoil sequential spe10

% The base case for the model is 2000 days. This can be reduced to make the
% problem faster to run.
[G, rock, fluid, deck, state] = setupEGG();
model = selectModelFromDeck(G, rock, fluid, deck);

schedule = convertDeckScheduleToMRST(model, deck);
% Set up the sequential model
seqModel = getSequentialModelFromFI(model, 'pressureLinearSolver', psolver,....
                                           'transportLinearSolver', tsolver);
% Run problem
[wsSeq, statesSeq, repSeq] = simulateScheduleAD(state, seqModel, schedule);

%% Simulate fully implicit
% Solve the fully implicit version of the problem, with a CPR
% preconditioner that uses the same linear solver for the pressure as the
% sequential solver.
solver = NonLinearSolver();
solver.LinearSolver = CPRSolverAD('ellipticSolver', psolver);

[wsFIMP, statesFIMP, repFIMP] = simulateScheduleAD(state, model, schedule, 'nonlinearsolver', solver);

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

%% Plot time taken
figure;
subplot(1, 2, 1)
bar([sum(repFIMP.SimulationTime), sum(repSeq.SimulationTime)]/60)
set(gca, 'XTickLabel', {'Fully-implicit', 'Sequential'})
ylabel('Simulation time [minutes]')
set(gca, 'FontSize', 12)
axis tight
subplot(1, 2, 2)
plotCellData(model.G, statesSeq{60}.s(:, 1), 'edgecolor', 'none')
title('S_w')
set(gca, 'FontSize', 12)

axis equal tight off
%% Copyright notice

% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>

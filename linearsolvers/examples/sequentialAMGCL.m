mrstModule add ad-core ad-blackoil spe10 sequential mrst-gui linearsolvers

% Select layer 1
layers = 1;

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

%%
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

function [states, failure] = runSimulationProblem(model, state0, schedule, varargin)
%Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    opt = struct('useCPR',  false, ...
                 'useAGMG', false, ...
                 'selectLinearSolver', false, ...
                 'stepcount',   inf, ...
                 'cprTol', 1e-3 ...
                );
   
    opt = merge_options(opt, varargin{:});
    
    if isfinite(opt.stepcount)
        schedule.step.val = schedule.step.val(1:opt.stepcount);
        schedule.step.control = schedule.step.control(1:opt.stepcount);
    end
    
    if opt.useCPR
        if opt.useAGMG
            mrstModule add agmg
            ellipSolver = AGMGSolverAD('tolerance', 1e-2');
        else
            ellipSolver = BackslashSolverAD();
        end
        linsolve = CPRSolverAD('ellipticSolver', ellipSolver, ...
                               'relativeTolerance', opt.cprTol);
    elseif opt.selectLinearSolver
        linsolve = selectLinearSolverAD(model);
    else
        linsolve = BackslashSolverAD();
    end
    % Set max substeps low because it is a failure if these tests do not
    % run through the entire schedule without cutting timesteps more than
    % once
    nonlinear = NonLinearSolver('LinearSolver', linsolve, ...
                                'maxTimestepCuts', 2,...
                                'errorOnFailure', false);
    % Ignore flux output for tests
    model.outputFluxes = false;
    
    [wellSols, states, report] = simulateScheduleAD(state0, model, schedule,...
                                'NonLinearSolver', nonlinear);
    failure = report.Failure;
end

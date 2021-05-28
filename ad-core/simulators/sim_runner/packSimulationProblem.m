function problem = packSimulationProblem(state0, model, schedule, BaseName, varargin)
%Pack simulation inputs into a single atomic representation of problem
%
% SYNOPSIS:
%   problem = packSimulationProblem(state0, model, schedule, 'MyCase')
%   problem = packSimulationProblem(state0, model, schedule, 'MyCase', 'Name', 'MyNewSolver', ...
%                                    'NonLinearSolver', nls)
%
% REQUIRED PARAMETERS:
%   initState    - Initial reservoir/model state. It should have whatever
%                  fields are associated with the physical model, with
%                  reasonable values. It is the responsibility of the user
%                  to ensure that the state is properly initialized.
%
%   model        - The physical model that determines jacobians/convergence
%                  for the problem. This must be a subclass of the
%                  `PhysicalModel` base class.
%
%   schedule     - Schedule containing fields step and control, defined as
%                  follows:
%                         - `schedule.control` is a struct array containing
%                           fields that the model knows how to process.
%                           Typically, this will be the fields such as `.W` 
%                           for wells or `.bc` for boundary conditions.
%
%                         - `schedule.step` contains two arrays of equal 
%                           size named `val` and `control`. Control is a
%                           index into the `schedule.control` array,
%                           indicating which control is to be used for the
%                           timestep.`schedule.step.val` is the timestep
%                           used for that control step.
%   BaseName      - Name of the case being simulated. This should be an
%                   unique identifier which can group several simulation
%                   realizations of the same case.
%
% OPTIONAL PARAMETERS:
%   NonLinearSolver - NonLinearSolver instance to be used for simulation.
%
%   Directory       - Folder where simulation output is to be stored. The
%                     default is fullfile(mrstOutputDirectory(), BaseName).
%
%   Name            - Name of this specific problem. We assume that
%                     BaseName + Name uniquely identifies all inputs to
%                     this function. Output will be stored under
%                     BaseName/Name.
%
%   Modules         - List of modules to be loaded before simulation.
%                     Defaults to the current loaded module set.
%
%   Description      - A description of the case. Can be useful to add
%                      additional information to the shorted Name +
%                      BaseName pair.
%
%   ExtraArguments   - Cell array of extra arguments to simulateScheduleAD.
%
% RETURNS:
%   problem - Self-contained packed simulation problem which can be passed
%             to a number of routines, e.g:
%                 - simulatePackedProblem: Simulate and store output.
%                 - getPackedSimulatorOutput: Retrieve results.
%                 - simulatePackedProblemBackground: Simulate in seperate
%                 thread.
% EXAMPLE:
%   demoPackedProblems
%
% SEE ALSO:
%   getPackedSimulatorOutput, simulatePackedProblem

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

    opt = struct('NonLinearSolver', [], ...
                 'Directory',          '', ...
                 'Name',            '', ...
                 'Modules',         {mrstModule()}, ...
                 'Description',     '', ...
                 'ExtraArguments', {{}} ...
                 );
    opt = merge_options(opt, varargin{:});
    noName = isempty(opt.Name);
    if noName
        opt.Name = class(model);
    end

    if isempty(opt.Description)
        if noName
            opt.Description = opt.Name;
        else
            opt.Description = [opt.Name, '_', class(model)];
        end
        if ~isempty(opt.NonLinearSolver)
            id = opt.NonLinearSolver.getId();
            if ~isempty(id)
                opt.Description = [opt.Description, '_', id];
            end
        end
    end

    if isempty(opt.Directory)
        opt.Directory = fullfile(mrstOutputDirectory(), BaseName);
    end
    
    problem.BaseName = BaseName;
    problem.Name = opt.Name;
    problem.Description = opt.Description;
    % OutputMinisteps needs to be treated specially by the packed problems.
    % We get the value and remove it from the extra arguments so it is
    % stored in a single location.
    simopt = struct('OutputMinisteps', false);
    [simopt, extra] = merge_options(simopt, opt.ExtraArguments{:});
    sim = struct('state0', state0, ...
                 'model', model, ...
                 'schedule', schedule, ...
                 'OutputMinisteps', simopt.OutputMinisteps, ...
                 'NonLinearSolver', opt.NonLinearSolver, ...
                 'ExtraArguments', {extra});
    problem.SimulatorSetup = sim;
    problem.Modules = opt.Modules;

    % We keep a separate prefix for substeps in case someone decides to
    % initialize a simulation problem twice with the same name and
    % different setting for output.
    if sim.OutputMinisteps
        simprefix = 'sub';
    else
        simprefix = '';
    end
    
    makeHandler = @(prefix) ResultHandler('dataPrefix', [simprefix, prefix], ...
                                          'writeToDisk', true,...
                                          'dataDirectory', opt.Directory, ...
                                          'dataFolder', problem.Name, ...
                                          'cleardir', false);
    

    problem.OutputHandlers.states   = makeHandler('state');
    problem.OutputHandlers.reports  = makeHandler('report');
    problem.OutputHandlers.wellSols = makeHandler('wellSols');
end
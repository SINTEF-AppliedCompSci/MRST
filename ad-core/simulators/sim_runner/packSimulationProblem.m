function problem = packSimulationProblem(state0, model, schedule, BaseName, varargin)
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
    
    sim = struct('state0', state0, ...
                 'model', model, ...
                 'schedule', schedule, ...
                 'NonLinearSolver', opt.NonLinearSolver, ...
                 'ExtraArguments', {opt.ExtraArguments});
    
    
    problem.SimulatorSetup = sim;
    problem.Modules = opt.Modules;

    
    makeHandler = @(prefix) ResultHandler('dataPrefix', prefix, ...
                                          'writeToDisk', true,...
                                          'dataDirectory', opt.Directory, ...
                                          'dataFolder', problem.Name, ...
                                          'cleardir', false);
    

    problem.OutputHandlers.states = makeHandler('state');
    problem.OutputHandlers.reports = makeHandler('report');
end
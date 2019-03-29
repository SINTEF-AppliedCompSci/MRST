function clearPackedSimulatorOutput(problems, varargin)
    opt = struct('Prompt', true,...
                 'Default', 'n', ...
                 'start', 1, ...
                 'stop', inf);
    opt = merge_options(opt, varargin{:});
    if isstruct(problems)
        problems = {problems};
    end
    for i = 1:numel(problems)
        problem = problems{i};
        s = problem.OutputHandlers.states;
        r = problem.OutputHandlers.reports;
        w = problem.OutputHandlers.wellSols;
        ns = s.numelData();
        nr = r.numelData();
        nw = w.numelData();
        if ns == 0 && nr == 0 && nw == 0
            continue
        end
        range = opt.start:min(opt.stop, ns);
        if opt.Prompt
            prompt = sprintf(['Do you want to delete %d states and ',...
                              'reports \nfor %s [%s]? y/n [%s]: '], ...
                              numel(range), problem.BaseName, problem.Name, opt.Default);
            str = input(prompt,'s');
            if ~strcmpi(str, 'y') && ~(strcmpi(opt.Default, 'y') && isempty(str))
                fprintf('Ok, will not remove files.\n');
                continue
            end
        end
        doPrint = mrstVerbose || opt.Prompt;
        dispif(doPrint, 'Removing files...');
        if opt.start == 1 && numel(range) == ns
            s.resetData();
            r.resetData();
            w.resetData();
        else
            s.resetData(range);
            r.resetData(range);
            w.resetData(range);
        end
        dispif(doPrint, ' Files removed.\n');
    end
end
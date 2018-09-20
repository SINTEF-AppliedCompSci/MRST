function clearPackedSimulatorOutput(problems, varargin)
    opt = struct('Prompt', true, 'Default', 'n');
    opt = merge_options(opt, varargin{:});
    if isstruct(problems)
        problems = {problems};
    end
    for i = 1:numel(problems)
        problem = problems{i};
        s = problem.OutputHandlers.states;
        r = problem.OutputHandlers.reports;
        
        ns = s.numelData();
        nr = r.numelData();
        if ns == 0 && nr == 0
            continue
        end
        if opt.Prompt
            prompt = sprintf(['Do you want to delete %d states and %d ',...
                              'reports \nfor %s [%s]? y/n [%s]: '], ...
                              ns, nr, problem.BaseName, problem.Name, opt.Default);
            str = input(prompt,'s');
            if ~strcmpi(str, 'y') && ~(strcmpi(opt.Default, 'y') && isempty(str))
                fprintf('Ok, will not remove files.\n');
                continue
            end
        end
        doPrint = mrstVerbose || opt.Prompt;
        dispif(doPrint, 'Removing files...');
        s.resetData();
        r.resetData();
        dispif(doPrint, ' Files removed.\n');
    end
end
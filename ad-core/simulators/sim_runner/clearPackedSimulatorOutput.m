function clearPackedSimulatorOutput(problems, varargin)
    opt = struct('Prompt', true);
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
                              'reports \nfor %s [%s]? Y/N [N]: '], ...
                              ns, nr, problem.BaseName, problem.Name);
            str = input(prompt,'s');
            if ~strcmpi(str, 'y')
                continue
            end
        end
        s.resetData();
        r.resetData();
    end
end
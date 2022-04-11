function [wells, states] = runJutulOnDaemon(state0, model, schedule, varargin)
    opt = struct('name', 'jutul_case', ...
                 'project', '', ...
                 'path', tempdir());
    [opt, extra] = merge_options(opt, varargin{:});
    if ~isempty(opt.project)
        opt.project = sprintf('--project="%s"', opt.project);
    end
    pth = writeJutulInput(state0, model, schedule, opt.name, 'printcmd', false, 'path', opt.path);
    % Create a julia script that runs the file
    %cmd_pth = fullfile(tempdir(), sprintf('run_%s.jl', opt.name));
    cmd_pth = fullfile(opt.path, sprintf('run_%s.jl', opt.name));
    fprintf('Creating Julia runscript at %s... ', cmd_pth) 
    f = fopen(cmd_pth, 'w');
    if f == 0
        error(ferror(f));
    end
    if mrstVerbose()
        info = 1;
    else
        info = -1;
    end
    fprintf('ok.\n');
    fprintf('Running Julia simulation...\n');
    % fprintf(f, 'println("Hello, world")');
    % fprintf(f, 'using LinearAlgebra\n');
    % fprintf(f, 'using Jutul\n');
    % fprintf(f, 'using JutulDarcy\n');
    fprintf(f, 'using JutulDarcy\nsimulate_mrst_case(\"%s\", write_mrst = true, info_level = %d, end_report = true, ascii_terminal = true)\n', ...
        pth, info);
    fclose(f);
    % Finally put together the command to invoke the daemon in client mode
    % and run the case.
    %cmd = sprintf('julia --startup-file=no --color=no -e "import Pkg; Pkg.add("DaemonMode"); using DaemonMode; runargs()" %s', cmd_pth);
    cmd = sprintf('julia %s --startup-file=no --color=no -e "using DaemonMode; runargs()" %s', ...
        opt.project, cmd_pth)
    id = system(cmd);
    if id == 1
        error('Something went wrong.')
    end
    fprintf('Julia simulation complete.\n');
    fprintf('Reading Julia output... ');
    [wells, states] = readJutulOutput(pth);
    fprintf('ok.\n');
end

function [wells, states] = runJutulOnDaemon(state0, model, schedule, varargin)
    opt = struct('name', 'jutul_case', ...
                 'project', '', ...
                 'path', tempdir());
    [opt, extra] = merge_options(opt, varargin{:});
    v = mrstVerbose();
    if ~isempty(opt.project)
        opt.project = sprintf('--project="%s"', opt.project);
    end
    pth = writeJutulInput(state0, model, schedule, opt.name, 'printcmd', false, 'path', opt.path);
    % Create a julia script that runs the file
    cmd_pth = fullfile(opt.path, sprintf('run_%s.jl', opt.name));
    dispif(v, 'Creating Julia runscript at %s... ', cmd_pth) 
    f = fopen(cmd_pth, 'w');
    if f == 0
        error(ferror(f));
    end
    if v
        info = 1;
    else
        info = -1;
    end
    extra_str = '';
    ne = numel(extra);
    assert(mod(ne, 2) == 0, 'Additional inputs must be key/value pairs')
    info_found = false;
    for i = 1:2:(ne-1)
        k = extra{i};
        if strcmpi(k, 'info_level')
            info_found = true;
        end
        val = extra{i+1};
        if isnumeric(val)
            val = num2str(val);
        elseif islogical(val)
            if val
                val = 'true';
            else
                val = 'false';
            end
        end
        extra_str = sprintf('%s, %s=%s', extra_str, k, val);
    end
    if ~info_found
        extra_str = sprintf('%s, info_level=%d', extra_str, info);
    end
    dispif(v, 'ok.\n');
    dispif(v, 'Running Julia simulation...\n');
    fprintf(f, 'using JutulDarcy\nsimulate_mrst_case(\"%s\", write_mrst = true%s, ascii_terminal = true)\n', ...
        pth, extra_str);
    fclose(f);
    % Finally put together the command to invoke the daemon in client mode
    % and run the case.
    cmd = sprintf('julia %s --startup-file=no --color=no -e "using DaemonMode; runargs()" %s', ...
        opt.project, cmd_pth);
    id = system(cmd);
    if id == 1
        error('Julia simulation was unable to complete successfully.')
    end
    dispif(v, 'Julia simulation complete.\n');
    dispif(v, 'Reading Julia output... ');
    [wells, states] = readJutulOutput(pth);
    dispif(v, 'ok.\n');
end

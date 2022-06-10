function [ws, states, reports] = simulateScheduleJutul(state0, model, schedule, varargin)
    opt = struct('name', class(model), 'daemon', true, 'pause', true, 'printcmd', true);
    opt = merge_options(opt, varargin{:});
    if opt.daemon
        [ws, states] = runJutulOnDaemon(state0, model, schedule, 'name', opt.name);
    else
        jpth = writeJutulInput(state0, model, schedule, opt.name, 'printcmd', opt.printcmd);
        if opt.pause
            disp('Pausing. Hit any key to continue once simulation has been run.')
            pause()
        end
        [ws, states] = readJutulOutput(jpth, 'error', true);
    end
    % Not parsed at the moment.
    reports = struct();
end

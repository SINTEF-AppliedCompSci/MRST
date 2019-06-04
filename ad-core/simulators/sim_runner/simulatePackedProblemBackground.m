function info = simulatePackedProblemBackground(problem, varargin)
    opt = struct('workdir', fullfile(mrstOutputDirectory(), 'sim_runner'), ...
                 'linux_arg', '-nodisplay', ...
                 'win_arg', '', ...
                 'matlab_arg', '-nosplash -noFigureWindows -nodesktop', ...
                 'extra_arg', '', ...
                 'verbose', mrstVerbose(), ...
                 'numCompThreads', -1);
    opt = merge_options(opt, varargin{:});

    basepath = fullfile(opt.workdir, [problem.BaseName, '_', problem.Name]);
    if not(exist(basepath, 'dir'))
        mkdir(basepath);
    end
    fn = makepath(basepath);
    logfile = makepath(fn, 'simulation.log');
    modlist = problem.Modules;
    save(makepath(fn, 'problem'), 'problem');
    save(makepath(fn, 'modlist'), 'modlist');
    matlab_options = build_matlab_options(opt);
    str = sprintf('matlab %s -sd "%s" -logfile %s -r "mrstModule add ad-core; simulatePackedProblemStandalone(''%s''); exit()" &', ...
        matlab_options, ...% Matlab options
        makepath(ROOTDIR()), ...% Startup dir - MRST root
        logfile, ... % Log file path
        fn ... % Problem path
        );
    disp(str);
    system(str);
    info = struct('path', basepath);
end

function options = build_matlab_options(opt)
    if ispc()
        % Windows arguments
        base = opt.win_arg;
    else
        % Mac OS and linux uses same parameters
        base = opt.linux_arg;
    end
    options = sprintf('%s %s %s', base, opt.matlab_arg, opt.extra_arg);
end

function str = makepath(varargin)
    str = fullfile(varargin{:});
    seperator = filesep();
    if strcmp(seperator, '\')
        str = strrep(str, '\', '\\');
    end
end

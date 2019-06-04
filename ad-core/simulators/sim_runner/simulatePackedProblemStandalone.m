function simulatePackedProblemStandalone(pth)
    p_path = fullfile(pth, 'problem.mat');
    % Load just the modules
    tmp = load(p_path, 'modlist');
    mrstModule('reset', tmp.modlist{:});
    % Now we can load the whole problem
    tmp = load(p_path, 'problem');
    lockpath = fullfile(pth, 'lock.mrst');
    if exist(lockpath, 'dir')
        error('This simulation is already running!')
    else
        fclose(fopen(lockpath, 'w'));
        try
            simulatePackedProblem(tmp.problem);
        catch

        end
        delete(lockpath);
    end
end

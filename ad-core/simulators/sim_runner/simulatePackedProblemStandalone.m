function simulatePackedProblemStandalone(pth)
    m_path = fullfile(pth, 'modlist.mat');
    tmp = load(m_path);
    mrstModule('reset', tmp.modlist{:});

    p_path = fullfile(pth, 'problem.mat');
    tmp = load(p_path);
    
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

function schedule = extend_migration_sch( schedule, new_mtime )
% Extend migration period of schedule

% Keeps same time step size. Only adds more steps to migration part
% (control=2) of schedule. Assuming schedule contains uniform migration
% step size.

% new_mtime is in seconds

    istepvec = schedule.step.val( schedule.step.control == 1 );
    isteps   = numel( schedule.step.control(schedule.step.control == 1) );
    mstepvec = schedule.step.val( schedule.step.control == 2 );
    
    
    % migration time
    tmp = cumsum(mstepvec);
    mtime = tmp(end); % s
    
    % number of migration steps
    msteps = numel(tmp);
    
    % migration time step size (uniform)
    dTm = mtime / msteps;
    

    % new number of migration steps (using dTm and new_mtime)
    new_msteps = new_mtime / dTm;
    
    % new mstepvec
    new_mstepvec = ones(new_msteps, 1) * dTm;
    
    % update schedule
    schedule.step.val       = [istepvec; new_mstepvec];
    schedule.step.control   = [ones(isteps, 1); ones(new_msteps, 1) * 2];
    

end


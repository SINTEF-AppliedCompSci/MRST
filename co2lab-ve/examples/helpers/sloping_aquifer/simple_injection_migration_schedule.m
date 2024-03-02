function schedule = simple_injection_migration_schedule(W, bc, inj_duration, ...
                                                        inj_steps, migration_duration, ...
                                                        migration_steps)

    % Define simple, single-well injection and migration scenario
    
    % Setting up two copies of the well and boundary specifications. 
    % Modifying the well in the second copy to have a zero flow rate.
    schedule.control    = struct('W', W, 'bc', bc);
    schedule.control(2) = struct('W', W, 'bc', bc);
    schedule.control(2).W.val = 0;

    dT_injection = rampupTimesteps(inj_duration, ...
                                   inj_duration/inj_steps, 7);
   
    dT_migration = repmat(migration_duration/migration_steps, migration_steps, 1);
    
    schedule.step.val = [dT_injection; dT_migration];
    schedule.step.control = [ones(numel(dT_injection), 1); ...
                             2 * ones(numel(dT_migration), 1)];

end

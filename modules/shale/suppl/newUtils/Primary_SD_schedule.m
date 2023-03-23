function schedule = Primary_SD_schedule(primary,SD, W)

    schedule = struct();
    [W_Primary, W_SD, W_equil] = deal(W);
    
    W_Primary(2).status = false;W_Primary(2).cstatus = false;
    W_equil(1).status = false;W_equil(2).status = false;
    W_equil(1).cstatus = false;W_equil(2).cstatus = false;
    
    schedule.control = [struct('W', W_Primary);...  % primary
                        struct('W', W_equil);
                        struct('W', W_SD)];     % SD;
    dt_primary = rampupTimesteps(primary(1), primary(2), primary(3)); 
    dt_equil = rampupTimesteps(40*day, 5*day, 7);
    dt_SD = rampupTimesteps(SD(1), SD(2), SD(3));
    
    dt = [];
    control_id = [1*ones(numel(dt_primary),1);2*ones(numel(dt_equil),1);3*ones(numel(dt_SD),1)];

    dt = [dt_primary;dt_equil;dt_SD];
    schedule.step.val = dt;
    schedule.step.control = control_id;
end


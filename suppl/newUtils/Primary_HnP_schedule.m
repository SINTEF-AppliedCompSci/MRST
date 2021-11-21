function schedule = Primary_HnP_schedule(primary,inj, soak, prod,runtime, W)
%Make a schedule for HnP with user-specified durations for injection,
%soaking, and production
%
% SYNOPSIS:
%   schedule = simpleSchedule(timesteps);
%   schedule = simpleSchedule(timesteps, 'W', W, 'src', src, 'bc', bc);
%
% PARAMETERS:
%   
%   inj     - array of [injection period(day) timestep size(day) rampup counts].
%   soak    - array of [soaking period(day) timestep size(day) rampup counts].
%   prod    - array of [production period(day) timestep size(day) rampup counts].
%   W -  Wells to be used in the schedule. The wells will be active in
%        all timesteps. W is a struct with the following order: injector,
%        soaking, and producer.
%
% RETURNS:
%   schedule - struct suitable for HnP 

    dt_cycle= inj(1) + soak(1) + prod(1);
    cycle_count = floor(runtime/dt_cycle);

    schedule = struct();
    [W_inj, W_soak, W_prod] = deal(W);
    
    W_inj(2).status = false;
    W_soak(1).status = false; W_soak(2).status = false;
    W_prod(1).status = false;
    
    schedule.control = [struct('W', W_inj);...  % injection
                        struct('W', W_soak);     % soaking
                        struct('W', W_prod)];... % production
    dt_primary = rampupTimesteps(primary(1), primary(2), primary(3)); 
    dt = [];
    control_id = [3*ones(numel(dt_primary),1)];
    dt_inj = rampupTimesteps(inj(1), inj(2), inj(3));
    dt_inj = [dt_inj(1);dt_inj(1);dt_inj];
    
    dt_soak = rampupTimesteps(soak(1), soak(2), soak(3));
    dt_prod = rampupTimesteps(prod(1), prod(2), prod(3));
    for i = 1:cycle_count
        dt = [dt;dt_inj;dt_soak;dt_prod];
        control_id = [control_id;[2;2;ones(numel(dt_inj(3:end)),1)];2*ones(numel(dt_soak),1);3*ones(numel(dt_prod),1)];
%         control_id = [control_id;ones(numel(dt_inj),1);2*ones(numel(dt_soak),1);3*ones(numel(dt_prod),1)];
    end
    dt = [dt_primary;dt];
    schedule.step.val = dt;
    schedule.step.control = control_id;
end


function dT = rampupTimesteps(time, dt, n)
    % time: total time
    % dt  : desired timesteps
    % n   : number of rampup steps
    if nargin < 3
        n = 8;
    end
    rampup = 2*dt;
    assert(time > rampup, 'Rampup time is larger than total time!');
    % Initial geometric series
    dt_init = (dt./2.^(n:-1:1))';
    % Remaining time that must be discretized
    dt_left = time - sum(dt_init);
    % Even steps
    dt_rem = repmat(dt, floor(dt_left/dt), 1);
    % Final ministep if present
    dt_final = time - sum(dt_init) - sum(dt_rem);
    if dt_final == 0
        dt_final = [];
    end
    % Combined timesteps
    dT = [dt_init; dt_rem; dt_final];
end
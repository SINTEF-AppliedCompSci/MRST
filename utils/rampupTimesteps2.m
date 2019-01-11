function dtVec = rampupTimesteps2(time, dt, n)

    if nargin < 3
        n = 8;
    end
    if time == 0
        dtVec = [];
        return
    end
    
    % Initial geometric series
    dtVec    = (dt./2.^[n n:-1:1])';
    isRampup = true(numel(dtVec),1);
    remTime  = time - sum(dtVec);
    % Fill in and/or crop/scale so that is fills up entire time
    while remTime ~= 0
        if remTime < 0
            ix = cumsum(dtVec) <= time;
            dtVec    = dtVec(ix);
            isRampup = isRampup(ix);
            remTime  = time - sum(dtVec);
            dtVec(isRampup) = dtVec(isRampup).*(sum(dtVec(isRampup)) + remTime)./sum(dtVec(isRampup));
        else
            nNew     = max(floor(remTime/dt),1);
            dtVec    = [dtVec; repmat(dt, nNew, 1)];
            isRampup = [isRampup; false(nNew,1)];
        end
        remTime  = time - sum(dtVec);
    end
    
%     dtFull = repmat(dt, floor(remTime/dt), 1);
%     
%     remTime = time - sum([dtRampup; dtFull]);
%     dtRampup = dtRampup.*(sum(dtRampup) + remTime)/sum(dtRampup);
%     dtVec = [dtRampup; dtFull];
% %     
%     cs_time = cumsum(dt_init);
%     if any(cs_time > time)
%         dt_init = dt_init(cs_time < time);
%     end
%     
%     % Remaining time that must be discretized
%     dt_left = time - sum(dt_init);
%     % Even steps
%     dt_rem = repmat(dt, floor(dt_left/dt), 1);
%     % Final ministep if present
%     dt_final = time - sum(dt_init) - sum(dt_rem);
%     % Less than to account for rounding errors leading to a very small
%     % negative time-step.
%     if dt_final <= 0
%         dt_final = [];
%     end
%     % Combined timesteps
%     dT = [dt_init; dt_rem; dt_final];
    
end
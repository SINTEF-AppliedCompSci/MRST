function [sts, maxPressDev_stepNum] = computeOverpressure(initState, sts, schedule)

    initPress   = initState.pressure;
    maxPressDev = 0;
    maxCO2sat   = 0;
    timeYr      = convertTo( cumsum(schedule.step.val), year);

    for i=1:numel(sts)

        sts{i}.s( sts{i}.s(:,2) < 0.00001, 2 ) =  nan;

        sts{i}.pressDev = sts{i}.pressure - initPress;

        % convert pressure deviation to bars
        sts{i}.pressDev = convertTo(sts{i}.pressDev, barsa);
        sts{i}.pressDev( sts{i}.pressDev < 0.5  ) = nan;

        % find max pressure dev and when it occurred
        maxPressDev_curr = max(sts{i}.pressDev);
        if maxPressDev_curr > maxPressDev
            maxPressDev          = maxPressDev_curr;
            maxPressDev_timeYr   = timeYr(i);
            maxPressDev_stepNum  = i;
        end

        % max co2 saturation and when it occurred
        maxCO2sat_curr = max(sts{i}.s(:,2));
        if maxCO2sat_curr > maxCO2sat
            maxCO2sat            = maxCO2sat_curr;
            maxCO2sat_timeYr     = timeYr(i);
            maxCO2sat_stepNum    = i;
        end
    end
    
    fprintf('\n Max press deviation of %d (bars) occurred at %d years (step %d) since sim start.\n', maxPressDev, maxPressDev_timeYr, maxPressDev_stepNum )
    fprintf('\n Max CO2 saturation of %d occurred at %d years (step %d) since sim start.\n', maxCO2sat, maxCO2sat_timeYr, maxCO2sat_stepNum )


%     for i = 1:numel(sts)
%        [m, l] =  max(sts{i}.pressDev);
%        maxOverpressure(i) = m;
%        locationOfMaxOP(i) = l;
%     end
%     [maxOP, locOfMaxOP] = max(maxOverpressure);
    
    
end


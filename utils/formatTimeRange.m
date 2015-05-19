function s = formatTimeRange(time)
% Small utility which returns a human readable string from seconds.
    s = '';

    timesys = {year, day, hour, second, second/1000};
    timen = {'Year', 'Day', 'Hour', 'Second', 'Millisecond'};
    added = false;
    for i = 1:numel(timesys)
        a = floor(time/timesys{i});
        if a > 0,
            if added
                space = ', ';
            else
                space = '';
            end

            plural = '';
            if a ~= 1, plural = 's'; end

            s = [s, sprintf([space, '%d ', timen{i}, plural], a)];     %#ok

            time  = time - a*timesys{i};
            added = true;
        end
    end
end

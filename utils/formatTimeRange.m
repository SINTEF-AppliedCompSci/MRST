function s = formatTimeRange(time, limit)
% Small utility which returns a human readable string from seconds.
    if nargin < 2
        % Second argument can be a limit of time units to output
        limit = inf;
    end
    s = '';

    timesys = {year, day, hour, second, second/1000};
    timen = {'Year', 'Day', 'Hour', 'Second', 'Millisecond'};
    added = false;
    count = 0;
    for i = 1:numel(timesys)
        a = floor(time/timesys{i});
        if a > 0,
            if added
                space = ', ';
            else
                space = '';
            end
            count = count + added;

            plural = '';
            if a ~= 1, plural = 's'; end
            if count == limit
                % Max depth reached, just print decimal
                s = [s, sprintf([space, '%1.2f ', timen{i}, plural], ...
                                                    time/timesys{i})]; %#ok
                break
            else
                s = [s, sprintf([space, '%d ', timen{i}, plural], a)]; %#ok
            end

            time  = time - a*timesys{i};
            added = true;
        end
    end
end

function s = formatTimeRange(time, limit)
% Small utility which returns a human readable string from seconds.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

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

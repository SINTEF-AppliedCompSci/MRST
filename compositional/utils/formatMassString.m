function s = formatMassString(mass, limit)
% Small utility which returns a human readable string from mass.

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
    if sign(mass) > 0
        s = '';
    elseif sign(mass) < 0
        s = '-';
    else
        s = '0 Kilogram';
    end
    mass = abs(mass);

    masssys = {mega*1000*kilogram, 1000*kilogram, kilogram, milli*kilogram};
    massn = {'Megatonne', 'Tonne', 'Kilogram', 'Gram'};
    added = false;
    count = 0;
    for i = 1:numel(masssys)
        a = floor(mass/masssys{i});
        if a > 0
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
                s = [s, sprintf([space, '%1.2f ', massn{i}, plural], ...
                                                    mass/masssys{i})]; %#ok
                break
            else
                s = [s, sprintf([space, '%d ', massn{i}, plural], a)]; %#ok
            end

            mass  = mass - a*masssys{i};
            added = true;
        end
    end
end

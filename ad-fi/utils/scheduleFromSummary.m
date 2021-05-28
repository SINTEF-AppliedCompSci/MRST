function schedule = scheduleFromSummary(schedule, summary)
    % There are n-1 control steps for n summary steps

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

    n_t = numel(schedule.step.control);
    schedule.control = repmat(schedule.control(1), n_t, 1);
    schedule.step.control = (1:n_t) .';

    for t = 1:n_t
        wconprod = schedule.control(t).WCONPROD;
        for w = 1:numel(wconprod(:, 1))
            name = wconprod(w,1);
            assert(strcmpi(wconprod(w, 3), 'bhp'));
            schedule.control(t).WCONPROD{w,9} = convertFrom(summary.get(name, 'WBHP', t+1), psia);
        end

        if 1
        wconinje = schedule.control(t).WCONINJE;
        for w = 1:numel(wconinje(:, 1))
            name = wconinje(w,1);
            assert(strcmpi(wconinje(w, 4), 'bhp'));

            schedule.control(t).WCONINJE{w,7} = convertFrom(summary.get(name, 'WBHP', t+1), psia);

            schedule.control(t).WCONINJE{w,5} = inf;
        end
        end
    end

end

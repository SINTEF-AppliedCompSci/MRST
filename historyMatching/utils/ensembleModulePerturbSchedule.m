function schedule = ensembleModulePerturbSchedule(schedule)
%Undocumented Utility Function

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    for n=2:max(schedule.step.control)
        schedule.control(n) = schedule.control(1);
        for i=1:numel(schedule.control(n).W)
            W = schedule.control(n).W(i);
            switch W.type
                case 'rate'
                    W.val = (.75 + .5*rand)*W.val;
                case 'bhp'
                    %if rand < 0.2
                    %    W.status = false;
                    % else
                        W.val = (.95 + 0.1*rand)*W.val;
                    %end
            end
            schedule.control(n).W(i) = W;
        end
    end
end

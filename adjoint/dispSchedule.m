function [] = dispSchedule(schedule, varargin)
%Display schedule values

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


opt = struct('Verbose' ,   mrstVerbose, ...
             'FileName',   [], ...
             'WellBased',  false);

opt = merge_options(opt, varargin{:});
wellBased = opt.WellBased;

numSteps = numel(schedule);

fprintf('\n----------------- DISPLAYING SCHEDULE ----------------\n')

if ~wellBased
    for k = 1 : numSteps
        fprintf('\nTime interval :  ');
        fprintf('%0.2f - %0.2f\n', schedule(k).timeInterval);
        fprintf('%9s%15s%15s\n', 'Well name', 'Type', 'Value')
        s = schedule(k);
        for k1 = 1:numel(s.names)
            fprintf('%9s%15s%15s\n', s.names{k1}, s.types{k1}, num2str(s.values(k1)))
        end
    end
else
    for k = 1 : numel(schedule(1).names)
        fprintf(['\nWell :  ', schedule(1).names{k}, '\n']);
        for k1 = 1:numSteps
            fprintf('%10.2f - %10.2f%10s%15s\n', schedule(k1).timeInterval, schedule(k1).types{k}, num2str(schedule(k1).values(k)));
        end
    end
end


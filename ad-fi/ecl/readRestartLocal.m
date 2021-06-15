function [rstrt, rsspec] = readRestartLocal(prefix, kwrds)
%Undocumented Utility Function

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

try
    rsspec = readEclipseOutputFileUnFmt([prefix, '.RSSPEC']);
catch me
    rsspec = [];
end

if (nargin < 2)
    if ~isempty(rsspec)
        kwrds = unique(rsspec.NAME.values);
        numSteps = numel(rsspec.TIME.values);
    else
        kwrds = {};
        numSteps = 1;
    end
else
    if (nargin == 2)
       numSteps = numel(rsspec.TIME.values);
    else

    end

end

fname = [prefix, '.UNRST'];

[fid, msg] = fopen(fname, 'r', 'ieee-be');
if fid < 0, error([fname, ': ', msg]); end

rstrt = struct();
for k = 1:numel(kwrds)
    rstrt.(genvarname(kwrds{k})) = repmat({{}}, [1, numSteps]);
end

fprintf(['Reading ', num2str(numel(kwrds)), ' fields in ', num2str(numSteps), ' steps\n'])
cnt = 0;
while ~feof(fid)
    [name, field] = readFieldUnFmt(fid);
    if strcmp(name, 'SEQNUM')
        fprintf('*');cnt = cnt+1;
    end
    if ismember(name, kwrds) || isempty(kwrds)
        vnm = genvarname(name);
        if ~isempty(field)
            rstrt.(vnm){cnt} = field.values;
        else
            rstrt.(vnm){cnt} = [];
        end
    end
end
fprintf('\n')
fclose(fid);



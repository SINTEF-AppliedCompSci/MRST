function printConvergenceReport(names, values, converged, iteration)
    nl = cellfun(@numel, names);
    sep = repmat('=', 1, sum(max(nl, 8)) + 3*numel(nl) + 8);
    if iteration == 1
        fprintf('%s\n', sep);
        fprintf('| It # ');
        fprintf('| %-8s ', names{:});
        fprintf('|\n%s\n', sep);
    end
    fprintf('| %4d ', iteration);
    for i = 1:numel(values)
        linen = max(nl(i), 8);
        fprintf(['| %-', num2str(linen), '.2e '], values(i));
    end
    fprintf('|\n')
    if converged
        fprintf('%s\n', sep);
    end
end

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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

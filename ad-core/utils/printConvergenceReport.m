function printConvergenceReport(names, values, converged, iteration, endOfBlock)
% Print a neatly formatted convergence report
%
% SYNOPSIS:
%   printConvergenceReport({'myEquation', 'yourEquation'}, [1, 25], [true, false], it);
%
% DESCRIPTION:
%   Print convergence report to the Command Window. Two lines are plotted
%   for the first iteration, and one line for succeeding iterations.
%
% REQUIRED PARAMETERS:
%   names     - Names of the different convergence measures. Cell array of
%               length N where N is the number of different measures (for
%               instance, residual norms for different equations)
%
%   values    - Double array of length N, where each entry corresponds to 
%               the current value of the different named measures.
%
%   converged - Boolean for each value indicating if convergence has been
%               achieved for that value.
%
% RETURNS:
%   Nothing.
%

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
    if nargin < 5
        endOfBlock = all(converged) && iteration > 1;
    end
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
        fprintf('|');
        if converged(i)
            fprintf('*');
        else
            fprintf(' ');
        end
        fprintf(['%-', num2str(linen), '.2e '], values(i));
    end
    fprintf('|\n')
    if all(endOfBlock)
        fprintf('%s\n', sep);
    end
end


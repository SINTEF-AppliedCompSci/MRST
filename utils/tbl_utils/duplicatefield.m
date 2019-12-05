function newtbl = duplicatefield(tbl, fdcell)
%
%
% SYNOPSIS:
%   function newtbl = duplicatefield(tbl, fdcell)
%
% DESCRIPTION: Duplicate an index
%
% PARAMETERS:
%   tbl    - Index table
%   fdcell - Names of the index to duplicate with names of the duplicated
%   indices. The syntax is {'name', {'dupname1', 'dupname2'}}
%
% RETURNS:
%   newtbl - The resulting index table
%
% SEE ALSO:
%   `setupTableMapping`

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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

    oldfd = fdcell{1};
    fd1 = fdcell{2}{1};
    fd2 = fdcell{2}{2};

    [a, fds] = convertTableToArray(tbl, {oldfd});
    a = [a(:, 1), a];
    fds = fds(2 : end);
    fds = {fd1, fd2, fds{:}};

    newtbl = convertArrayToTable(a, fds);
end

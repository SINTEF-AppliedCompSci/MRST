function projtbl = projTable(tbl, fds)
%
%
% SYNOPSIS:
%   function projtbl = projTable(tbl, fds)
%
% DESCRIPTION: Project an index table on a subset of indices, meaning that
% given a table tbl with field tbl.A and tbl.B. If we project on the index
% name 'A', we construct a projtbl with field projtbl.A such that, for all i,
% there exists a unique j such that tbl.A(i) = projtbl.A(j) 
%
% PARAMETERS:
%   tbl - Index table
%   fds - field names along which we project. Syntax is {'A1', 'A2'} for example
%
% RETURNS:
%   projtbl - resulting index table
%
% SEE ALSO:
%   `setupTableMapping`.

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

    a = convertTableToArray(tbl, fds);
    nfds = numel(fds);
    a = a(:, (1 : nfds));
    a = unique(a, 'rows');
    projtbl = convertArrayToTable(a, fds);
end

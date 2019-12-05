function tbl = convertArrayToTable(A, fds)
%
%
% SYNOPSIS:
%   function tbl = convertArrayToTable(A, fds)
%
% DESCRIPTION: Convert an array to an indexing table. The indexing table are
% given as structure with field names. Each field contains a vector of
% indices. All these vectors have same length. 
%
% PARAMETERS:
%   A   - Array
%   fds - Field names given to each of the column of A
%
% RETURNS:
%   tbl - The corresponding indexing table
%
% SEE ALSO:
%   `convertTableToArray`, `setupTableMapping`.

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

    sz = size(A);
    nfds = numel(fds);
    
    for ifield = 1 : nfds
        tbl.(fds{ifield}) = A(:, ifield);
    end
    tbl.num = sz(1);
end

function newtbl = replacefield(tbl, fieldpairs)
%
%
% SYNOPSIS:
%   function newtbl = replacefield(tbl, fieldpairs)
%
% DESCRIPTION: Utility function to replace field names in an index table
%
% PARAMETERS:
%   tbl        - Index table
%   fieldpairs - The syntax is either fieldpairs = {'A', 'B'} to replace name
%   'A' with 'B' or fieldpairs = {{'A1', 'B1'}, {'A2', 'B2'}} to replace
%   names 'A1' and 'A2' with 'B1' and 'B2', respectively
%
% RETURNS:
%   newtbl - Resulting index table
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

    newtbl = tbl;
    if iscell(fieldpairs{1})
        for i = 1 : numel(fieldpairs)
            newtbl = replacethisfield(newtbl, fieldpairs{i});
        end
    else
        fieldpair = fieldpairs;
        newtbl = replacethisfield(newtbl, fieldpair);
    end
end

function tbl = replacethisfield(tbl, fieldpair)
    oldfield = fieldpair{1};
    newfield = fieldpair{2};
    tbl.(newfield) = tbl.(oldfield);
    tbl = rmfield(tbl, oldfield);
end

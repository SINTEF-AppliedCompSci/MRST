function tbl = addLocInd(tbl, locindname)
%
%
% SYNOPSIS:
%   function tbl = addLocInd(tbl, locindname)
%
% DESCRIPTION: Add a new field with a local indexing (1,2, ..., size of
% table) to the given indexing table.
%
% PARAMETERS:
%   tbl        - Indexing table
%   locindname - Chosen name for the local index
%
% RETURNS:
%   tbl - Indexing table with new field giving a local index (from 1 to size
%   of table)

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

    tbl.(locindname) = (1 : tbl.num)';
end

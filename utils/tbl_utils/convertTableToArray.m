function [A, fds] = convertTableToArray(tbl, fdsfirst)
%
%
% SYNOPSIS:
%   function [A, fds] = convertTableToArray(tbl, fdsfirst)
%
% DESCRIPTION: Convert an indexing table to an array. The indexing table are
% given as structure with field names. Each field contains a vector of
% indices. All these vectors have same length.
%
% PARAMETERS:
%   tbl      - Indexing table
%   fdsfirst - We give the fields that we want to appear first and in the
%   given order in the array
%
% RETURNS:
%   A   - Corresponding array
%   fds - returns the name of the fields corresponding to the columns in A
%
% SEE ALSO:
%   `convertArrayToTable`, `setupTableMapping`.

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
    
    fds = fieldnames(tbl);
    
    fdsToRemove = {fdsfirst{:}, 'num'};
    for ifield = 1 : numel(fdsToRemove)
        fds = fds(~strcmp(fds, fdsToRemove{ifield}));
    end
    
    fds = {fdsfirst{:}, fds{:}};
    A = [];
    A = tbl.(fds{1});
    for ifield = 2 : numel(fds)
        A = [A, tbl.(fds{ifield})];
    end
     
end

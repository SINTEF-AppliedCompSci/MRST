function tbl = sortTable(tbl, fds, varargin)
%
%
% SYNOPSIS:
%   tbl = sortTable(tbl, fds, varargin)
%
% DESCRIPTION:
%   Sort in increasing order the rows of indices given by the indexing
%   table tbl, following the order given by fds. (We use sortrows)
%
% PARAMETERS:
%   tbl      - Indexing table
%   fds      - Field names which give the order of the columns for the row
%              indexing sort
%   varargin - 
%
% KEYWORD ARGUMENTS:
%   keepAllFields - We keep all the fields for the table, even those that are
%                   not given in fds but belong to tbl. The default value is false. This
%                   option has not been properly tested and we recommend not using it for the moment.
%
% RETURNS:
%   tbl - Indexing table that has been sorted by row.
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

    opt = struct('keepAllFields', false);
    opt = merge_options(opt, varargin{:});
    
    a = convertTableToArray(tbl, fds);
    a = sortrows(a);
    
    if opt.keepAllFields
        warning('option not really tested yet!');
        ofds = fieldnames(tbl);
    
        fdsToRemove = {fds{:}, 'num'};
        for ifield = 1 : numel(fdsToRemove)
            ofds = ofds(~strcmp(ofds, fdsToRemove{ifield}));
        end
        
        fds = {fds{:}, ofds{:}};
    end
    
    tbl = convertArrayToTable(a, fds);
end

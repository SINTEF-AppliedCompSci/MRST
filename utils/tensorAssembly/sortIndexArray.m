function tbl = sortIndexArray(tbl, fds, varargin)
%
%
% SYNOPSIS:
%   tbl = sortIndexArray(tbl, fds, varargin)
%
% DESCRIPTION: . We use function sortrows to sort in increasing order the rows
% of indices given by the IndexArray tbl, following the order given by fds.
%
% PARAMETERS:
%   tbl      - IndexArray
%   fds      - Field names which give the order of the columns for the row
%              indexing sort
%   varargin - 
%
% KEYWORD ARGUMENTS:
%   keepAllFields - We keep all the fields for the table, even those that are
%                   not given in fds but belong to tbl. The default value is
%                   false. This option has not been properly tested and we
%                   recommend not using it for the moment.
%
% RETURNS:
%   tbl - IndexArray that has been sorted by rows.
%
% SEE ALSO: `IndexArray`, `crossIndexArray`, `projIndexArray`
%   `crossIndexArray`, `IndexArray`.

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

    opt = struct('keepAllFields', false);
    opt = merge_options(opt, varargin{:});
    
    fdnames = tbl.fdnames;
    [isfd, fdind2] = ismember(fds, fdnames);
    
    assert(all(isfd), 'fields not found in IndexArray');
    
    inds = tbl.inds;
    inds = inds(:, fdind2);
    
    inds = sortrows(inds);
    fdnames = fds;
    
    
    if opt.keepAllFields
        warning('option not really tested yet!');
        ofds = fieldnames(tbl);
    
        fdsToRemove = {fds{:}, 'num'};
        for ifield = 1 : numel(fdsToRemove)
            ofds = ofds(~strcmp(ofds, fdsToRemove{ifield}));
        end
        
        fds = {fds{:}, ofds{:}};
    end

    tbl.fdnames = fdnames;
    tbl.inds = inds;
    
end

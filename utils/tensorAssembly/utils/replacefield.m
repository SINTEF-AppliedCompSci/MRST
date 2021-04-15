function newtbl = replacefield(tbl, fieldpairs)
%
%
% SYNOPSIS:
%   function newtbl = replacefield(tbl, fieldpairs)
%
% DESCRIPTION: Utility function to replace field names in an IndexArray
%
% PARAMETERS: tbl - IndexArray fieldpairs - The syntax is either fieldpairs =
%   {'A', 'B'} to replace name 'A' with 'B' or fieldpairs = {{'A1', 'B1'},
%   {'A2', 'B2'}} to replace names 'A1' and 'A2' with 'B1' and 'B2',
%   respectively. The new names (here 'B', 'B1' or 'B2') cannot coincide with
%   existing fieldnames of the IndexArray.
%
%   If a name pair is given the extra argument 'interchange', as for example
%   {'A', 'B', 'interchange'}, then the effect is to interchange the fields
%   'A' and 'B'.
%  
%   If the second argument is the empty string, for example {'A', ''}, then
%   the field is removed
%
%   newtbl - Resulting IndexArray
%
% SEE ALSO:
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
    
    oldfdname = fieldpair{1};
    newfdname = fieldpair{2};
    interchange = false;
    if (numel(fieldpair) > 2)
        assert(strcmp('interchange', fieldpair{3}), 'do not recognized third argument');
        interchange = true;
    end
    
    fdnames = tbl.fdnames;
    n = numel(fdnames);
    
    [isold, oldfdind] = ismember(oldfdname, fdnames);
    assert(isold, 'fieldname not recognized in IndexArray');
    
    if interchange
        [isnew, newfdind] = ismember(newfdname, fdnames);
        assert(isnew, 'fieldname not recognized in IndexArray');
        fdnames{oldfdind} = newfdname;
        fdnames{newfdind} = oldfdname;
    else
        if ~strcmp(newfdname, '')
            fdnames{oldfdind} = newfdname;
        else
            % we erase this field
            fdinds = (1 : n)';
            fdinds(oldfdind) = [];
            fdnames = fdnames(fdinds);
            tbl.inds = tbl.inds(:, fdinds);
        end            
    end
    
    tbl.fdnames = fdnames;
    
end

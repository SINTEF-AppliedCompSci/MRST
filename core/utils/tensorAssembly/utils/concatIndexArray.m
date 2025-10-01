function tbl = concatIndexArray(tbl1, tbl2, fdnames, varargin)
%Undocumented Utility Function

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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
    
    opt = struct('checkUnique', true, ...
                 'removeDuplicates', false);
    opt = merge_options(opt, varargin{:}); 

    if nargin < 3 || isempty(fdnames)
        fdnames1 = tbl1.fdnames;
        fdnames2 = tbl2.fdnames;
        assert(all(ismember(fdnames1, fdnames2)) && all(ismember(fdnames1, fdnames2)), 'not matching field names');
        fdnames = fdnames1;
    end

    if tbl1.isvirtual || tbl2.isvirtual
        assert(tbl1.isvirtual && tbl2.isvirtual, 'in case of virtual array, we allow only for concatenation if both indexarrays are virtual');
        tbl = IndexArray([], 'fdnames', fdnames, 'isvirtual', true, 'num', tbl1.num + tbl2.num);
        if opt.checkUnique
            warning('cannot check uniqueness of indices in case of virtual index arrays');
        end
        if opt.removeDuplicates
            warning('we cannot remove duplicates in case of virtual index arrays');
        end
        
        return
        
    end
    
    inds1 = tbl1.gets(fdnames);
    inds2 = tbl2.gets(fdnames);

    inds = [inds1; inds2];
    
    if (opt.checkUnique) & ~opt.removeDuplicates
        indstest = unique(inds, 'rows');
        assert(size(indstest, 1) == size(inds, 1), 'there are repeated indices');
    end    
    
    if opt.removeDuplicates
        inds = unique(inds, 'rows');
    end
    
    tbl = IndexArray([], 'inds', inds, 'fdnames', fdnames);
end

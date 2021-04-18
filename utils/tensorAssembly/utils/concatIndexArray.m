function tbl = concatIndexArray(tbl1, tbl2, fdnames, varargin)
%Undocumented Utility Function

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
    
    opt = struct('checkUnique', true, ...
                 'removeDuplicates', false);
    opt = merge_options(opt, varargin{:}); 
    
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

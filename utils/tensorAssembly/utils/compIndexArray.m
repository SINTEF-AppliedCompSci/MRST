function iseq = compIndexArray(tbl1, tbl2)
%
%
% SYNOPSIS:
%   iseq = compIndexArray(tbl1, tbl2)
%
% DESCRIPTION:
%   Compares the two IndexArrays tbl1 and tbl2 (field names and indices)
%
% PARAMETERS:
%   tbl1 - IndexArray
%   tbl2 - IndexArray
%
% RETURNS:
%   iseq - true if equal false otherwise

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

    fds1 = tbl1.fdnames;
    fds2 = tbl2.fdnames;    

    if tbl1.isvirtual || tbl2.isvirtual
        % the comparison cannot be done since the Index arrays are virtual
        iseq = NaN;
        return
    end

    iseq = (numel(fds1) == numel(fds2));

    if iseq

        [~, fdind2] = ismember(fds1, fds2);

        mat1 = sortrows(tbl1.inds);
        mat2 = sortrows(tbl1.inds(:, fdind2));

        iseq = all(all((mat1 - mat2) == 0));
    end
end

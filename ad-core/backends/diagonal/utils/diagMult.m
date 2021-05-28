function [x, D] = diagMult(v, M, D)
%Internal function for diagonal multiplication in AD code

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

    if ~any(v)
        sz = size(M);
        x = sparse([],[],[],sz(1), sz(2));
    elseif nnz(M) == 0
        x = M;
    elseif isscalar(v)
        x = v*M;
    else
        if isempty(D)
            n = numel(v);
            ix = (1:n)';
            D = sparse(ix, ix, v, n, n);
        end
        x = D*M;
    end
end

function Tc = upscaleByMsBasis(CG, T, varargin)
%Undocumented Utility Function

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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

    G = CG.parent;
    A = getIncomp1PhMatrix(G, T);
    basis = getMultiscaleBasis(CG, A, varargin{:});
    
    A_c = basis.R*A*basis.B;
    
    [ic, jc, vc] = find(A_c);
    Tc = ones(CG.faces.num, 1)*sqrt(eps)*max(abs(vc));

    for i = 1:CG.cells.num
        fa = gridCellFaces(CG, i);
        for j = 1:numel(fa)
            f = fa(j);
            
            c = CG.faces.neighbors(f, :);
            c = c(c~=i);
            if c == 0
                continue
            end
            Tc(f) = abs(vc(ic == i & jc == c));
        end
    end
end

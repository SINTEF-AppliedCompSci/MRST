function v = singlePointUpwind(flag, N, v, useMex)
%Single-point upwind for the GenericAD library

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

    vD = value(v);
    if useMex
        val = mexSinglePointUpwindVal(vD, N, flag);
    else
        cells = N(:, 2);
        cells(flag) = N(flag, 1);
        val = vD(cells, :);
    end
    if isa(v, 'GenericAD')
        M = [];
        DS = [];
        v.val = val;
        for jacNo = 1:numel(v.jac)
            [v.jac{jacNo}, M, DS] = upwindJac(v.jac{jacNo}, flag, N, M, DS, useMex);
        end
    else
        v = val;
    end
end

function [jac, M, DS] = upwindJac(jac, flag, N, M, DS, useMex)
    if issparse(jac)
        if any(jac(:))
            if isempty(M)
                nf = size(N, 1);
                nc = max(max(N));
                upCell = flag.*N(:, 1) + ~flag.*N(:, 2);
                M = sparse((1:nf)', upCell, 1, nf, nc);
            end
            jac = M*jac;
        else
            jac = sparse([], [], [], size(N, 1), matrixDims(jac, 2));
        end
    else
        if jac.isZero
            jac = jac.toZero(size(N, 1));
        else
            if useMex
                diagonal = mexSinglePointUpwindDiagonalJac(jac.diagonal, N, flag);
            else
                nf = size(N, 1);
                nder = size(jac.diagonal, 2);
                diagonal = zeros(nf, 2*nder);
                diagonal(flag, 1:nder) = jac.diagonal(N(flag, 1), :);

                notFlag = ~flag;
                diagonal(notFlag, (nder+1):end) = jac.diagonal(N(notFlag, 2), :);
                jac.diagonal = diagonal;
            end

            if isempty(DS)
                jac = DiagonalSubset(diagonal, jac.dim, N, [], jac.subset);
                DS = jac;
            else
                DS.diagonal = diagonal;
                DS.map = N;
                DS.dim = jac.dim;
                DS.parentSubset = jac.subset;
                jac = DS;
            end
        end
    end

end

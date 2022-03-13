function v = twoPointGradient(N, v, M, useMex)
%Discrete gradient for the GenericAD library

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    if isa(v, 'GenericAD')
        nc = numel(v.val);
        v.val = gradVal(v.val, N, useMex);
        v.jac = cellfun(@(x) gradJac(x, N, nc, useMex), v.jac, 'UniformOutput', false);
    else
        v = gradVal(v, N, useMex);
    end
end

function jac = gradJac(jac, N, nc, useMex)
    if issparse(jac)
        if nnz(jac)
            jac = jac(N(:, 2), :) - jac(N(:, 1), :);
        else
            jac = sparse([], [], [], size(N, 1), size(jac, 2));
        end
    elseif jac.isZero
        jac = jac.toZero(size(N, 1));
        return
    else
        rowMajor = jac.rowMajor;
        if useMex
            diagonal = mexTwoPointGradientDiagonalJac(jac.diagonal, N, nc, rowMajor);
        else
            if rowMajor
                diagonal = [-jac.diagonal(:, N(:, 1)); jac.diagonal(:, N(:, 2))];
            else
                diagonal = [-jac.diagonal(N(:, 1), :), jac.diagonal(N(:, 2), :)];
            end
        end
        jac = FixedWidthJacobian(diagonal, jac.dim, N, [], jac.subset, useMex, rowMajor, 'interiorfaces');
    end
end

function v = gradVal(val, N, useMex)
    if useMex
        v = mexTwoPointGradientVal(val, N);
    else
        if size(val, 2) > 1
            v = val(N(:, 2), :) - val(N(:, 1), :);
        else
            v = diff(val(N), 1, 2);
        end
    end
end

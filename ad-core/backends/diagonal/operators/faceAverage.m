function v = faceAverage(N, v, useMex)
% Face average operator for the GenericAD library

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
    if isa(v, 'GenericAD')
        nc = numel(v.val);
        if useMex
            v.val = mexFaceAverageVal(v.val, N);
        else
            if size(N, 1) == 1
                v.val = 0.5*sum(v.val(N), 1);
            else
                v.val = 0.5*sum(v.val(N), 2);
            end
        end
        v.jac = cellfun(@(x) avgJac(x, N, nc, useMex), v.jac, 'UniformOutput', false);
    elseif useMex
        v = mexFaceAverageVal(v, N);
    else
        v = 0.5*(v(N(:, 1), :) + v(N(:, 2), :));
    end
end

function jac = avgJac(jac, N, nc, useMex)
    if issparse(jac)
        if any(jac(:))
            jac = 0.5*(jac(N(:, 1), :) + jac(N(:, 2), :));
        else
            jac = sparse([], [], [], size(N, 1), size(jac, 2));
        end
    elseif jac.isZero
        jac = jac.toZero(size(N, 1));
    else
        rowMajor = jac.rowMajor;
        if useMex
            diagonal = mexFaceAverageDiagonalJac(jac.diagonal, N, nc, rowMajor);
        else
            if rowMajor
                diagonal = 0.5*[jac.diagonal(:, N(:, 1)); jac.diagonal(:, N(:, 2))];
            else
                diagonal = 0.5*[jac.diagonal(N(:, 1), :), jac.diagonal(N(:, 2), :)];
            end
        end
        jac = FixedWidthJacobian(diagonal, jac.dim, N, [], jac.subset, useMex, rowMajor, 'interiorfaces');
    end
end

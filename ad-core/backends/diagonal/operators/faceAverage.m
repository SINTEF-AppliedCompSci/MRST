function v = faceAverage(N, v, useMex)
% Face average operator for the GenericAD library

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

    if isa(v, 'GenericAD')
        if size(N, 1) == 1
            v.val = 0.5*sum(v.val(N), 1);
        else
            v.val = 0.5*sum(v.val(N), 2);
        end
        v.jac = cellfun(@(x) avgJac(x, N, useMex), v.jac, 'UniformOutput', false);
    else
        v = 0.5*(v(N(:, 1), :) + v(N(:, 2), :));
    end
end

function jac = avgJac(jac, N, useMex)
    if issparse(jac)
        if any(jac(:))
            jac = 0.5*(jac(N(:, 1), :) + jac(N(:, 2), :));
        else
            jac = sparse([], [], [], size(N, 1), size(jac, 2));
        end
    elseif jac.isZero
        jac = jac.toZero(size(N, 1));
    else
        if useMex
            diagonal = mexFaceAverageDiagonalJac(jac.diagonals, N);
        else
            diagonal = 0.5*jac.diagonals(:, N);
        end
        jac = DiagonalSubset(diagonal, jac.dim, N, [], jac.subset);
    end

end

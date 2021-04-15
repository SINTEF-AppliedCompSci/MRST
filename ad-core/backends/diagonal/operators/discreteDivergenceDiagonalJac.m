function jac = discreteDivergenceDiagonalJac(acc, jac, opt, getConservation)
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

    if nargin < 4
        getConservation = opt.useConservationJac;
    end
    
    if getConservation
        jac = ConservationLawJacobian(acc, jac, opt);
    else
        if isempty(acc)
            jac = divJac(jac, opt);
        else
            jac = accDivJac(acc, jac, opt);
        end
    end
end

function jac = divJac(jac, opt)
    if issparse(jac)
        if nnz(jac) > 0
            if isempty(opt.C)
                nf = opt.nf;
                nc = opt.nc;

                opt.C  = sparse(opt.N, [(1:nf)'; (1:nf)'], ones(nf,1)*[1 -1], nc, nf);
            end
            jac = opt.C*jac;
        else
            jac = sparse([], [], [], opt.nc, matrixDims(jac, 2));
        end
    elseif jac.isZero
        jac = sparse([], [], [], opt.nc, prod(jac.dim));
        return
    else
        if opt.useMex && (isempty(jac.parentSubset) || all(jac.parentSubset == (1:jac.dim(1))'))
            p = opt.mex;
            jac = mexDiscreteDivergenceJac([], jac.diagonal, opt.N, p.facePos, p.faces, p.cells, p.cellIndex, jac.rowMajor);
        else
            jac = opt.sortIx.C*jac.sparse();
        end
    end
end

function jac = accDivJac(acc, jac, opt)
    if issparse(jac)
        if nnz(jac) > 0
            if isempty(opt.C)
                nf = opt.nf;
                nc = opt.nc;
                opt.C = sparse(N, [(1:nf)'; (1:nf)'], ones(nf,1)*[1 -1], nc, nf);
            end
            jac = opt.C*jac + acc;
        else
            jac = acc;
        end
    elseif jac.isZero
            jac = acc;
        return
    else
        if opt.useMex && (isempty(jac.parentSubset) || (numel(jac.parentSubset) == jac.dim(1)) && all(jac.parentSubset == (1:jac.dim(1))'))
            p = opt.mex;
            if isa(acc, 'DiagonalJacobian')
                % NB currently not checking subset here - bug
                jac = mexDiscreteDivergenceJac(acc.diagonal, jac.diagonal, opt.N, p.facePos, p.faces, p.cells, p.cellIndex, jac.rowMajor);
            else
                jac = acc + mexDiscreteDivergenceJac([], jac.diagonal, opt.N, p.facePos, p.faces, p.cells, p.cellIndex, jac.rowMajor);
            end
        else
            jac = acc + opt.sortIx.C*jac.sparse();
        end
    end
end

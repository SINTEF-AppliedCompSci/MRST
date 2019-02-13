function v = discreteDivergence(N, v, nc, nf, sortIx, C, prelim, useMex)
% Discrete divergence for the NewAD library
    if isa(v, 'NewAD')
        for i = 1:numel(v.jac)
            v.jac{i} = divJac(v.jac{i}, N, nc, nf, sortIx, C, prelim, useMex);
        end
        v.val = accumulate(N, value(v), nc);
    else
        v = accumulate(N, value(v), nc);
    end
end

function v = accumulate(N, v, nc)
    v = accumarray(N(:, 1), v, [nc, 1]) - accumarray(N(:, 2), v, [nc, 1]);
end

function jac = divJac(jac, N, nc, nf, sortIx, C, prelim, useMex)
    useDivTerm = false;
    if issparse(jac)
        if nnz(jac) > 0
            if isempty(C)
                C  = sparse(N, [(1:nf)'; (1:nf)'], ones(nf,1)*[1 -1], nc, nf);
            end
            jac = C*jac;
        else
            if useDivTerm
                jac = DivergenceTerm([], [], [], nc, prod(jac.dim));
            else
                jac = sparse([], [], [], nc, matrixDims(jac, 2));
            end
        end
    elseif jac.isZero
%         jac = DiagonalJacobian();
%         if useDivTerm
%             jac = DivergenceTerm([], [], [], nc, prod(jac.dim));
%         else
            jac = sparse([], [], [], nc, prod(jac.dim));
%         end
        return
    else
        if useDivTerm
            jac = DivergenceTerm(jac, N, sortIx.C, prelim, useMex);
        else
            % Sparse version
            if useMex && (isempty(jac.parentSubset) || all(jac.parentSubset == (1:jac.dim(1))'))
                jac = mexDiscreteDivergenceJac(jac.diagonal, N, prelim.facePos, prelim.faces, prelim.cells, prelim.cellIndex);
            else
                jac = sortIx.C*jac.sparse();
            end
        end
    end
end
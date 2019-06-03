function v = discreteDivergence(acc, N, v, nc, nf, sortIx, C, prelim, useMex)
% Discrete divergence for the GenericAD library
    hasAcc = not(isempty(acc));
    if isa(v, 'GenericAD')
        v.val = accumulate(N, value(v), nc);
        if hasAcc && isa(acc, 'GenericAD')
            % Both present, both are AD
            for i = 1:numel(v.jac)
                v.jac{i} = accDivJac(acc.jac{i}, v.jac{i}, N, nc, nf, sortIx, C, prelim, useMex);
            end
            v.val = v.val + acc.val;
        else
            for i = 1:numel(v.jac)
                v.jac{i} = divJac(v.jac{i}, N, nc, nf, sortIx, C, prelim, useMex);
            end
            if hasAcc
                % Acc is not AD
                v.val = v.val + acc;
            end
        end
    else
        assert(isnumeric(v), 'Expected numeric vector, but got ''%s''\n', class(v))
        v = accumulate(N, v, nc);
        if hasAcc
            v = v + acc;
        end
    end
end

function v = accumulate(N, v, nc)
    v = accumarray(N(:, 1), v, [nc, 1]) - accumarray(N(:, 2), v, [nc, 1]);
end

function jac = divJac(jac, N, nc, nf, sortIx, C, prelim, useMex)
    if issparse(jac)
        if nnz(jac) > 0
            if isempty(C)
                C  = sparse(N, [(1:nf)'; (1:nf)'], ones(nf,1)*[1 -1], nc, nf);
            end
            jac = C*jac;
        else
            jac = sparse([], [], [], nc, matrixDims(jac, 2));
        end
    elseif jac.isZero
            jac = sparse([], [], [], nc, prod(jac.dim));
        return
    else
        if useMex && (isempty(jac.parentSubset) || all(jac.parentSubset == (1:jac.dim(1))'))
            jac = mexDiscreteDivergenceJac([], jac.diagonal, N, prelim.facePos, prelim.faces, prelim.cells, prelim.cellIndex);
        else
            jac = sortIx.C*jac.sparse();
        end
    end
end

function jac = accDivJac(acc, jac, N, nc, nf, sortIx, C, prelim, useMex)
    if issparse(jac)
        if nnz(jac) > 0
            if isempty(C)
                C  = sparse(N, [(1:nf)'; (1:nf)'], ones(nf,1)*[1 -1], nc, nf);
            end
            jac = C*jac + acc;
        else
            jac = acc;
        end
    elseif jac.isZero
            jac = acc;
        return
    else
        if useMex && (isempty(jac.parentSubset) || all(jac.parentSubset == (1:jac.dim(1))'))
            if isa(acc, 'DiagonalJacobian')
                % NB currently not checking subset here - bug
                jac = mexDiscreteDivergenceJac(acc.diagonal, jac.diagonal, N, prelim.facePos, prelim.faces, prelim.cells, prelim.cellIndex);
            else
                jac = acc + mexDiscreteDivergenceJac([], jac.diagonal, N, prelim.facePos, prelim.faces, prelim.cells, prelim.cellIndex);
            end
        else
            jac = acc + sortIx.C*jac.sparse();
        end
    end
end
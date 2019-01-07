function v = singlePointUpwind(flag, N, v, useMex)
    % Single-point upwind for the NewAD library
    vD = double(v);
    if useMex
        value = mexSinglePointUpwindVal(vD, N, flag);
    else
        cells = N(:, 2);
        cells(flag) = N(flag, 1);
        value = vD(cells, :);
    end
    if isa(v, 'NewAD')
        M = [];
        DS = [];
        v.val = value;
        for jacNo = 1:numel(v.jac)
            [v.jac{jacNo}, M, DS] = upwindJac(v.jac{jacNo}, flag, N, M, DS, useMex);
        end
    else
        v = value;
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
            elseif 1
                nf = size(N, 1);
                diagonal = zeros(2*nf, size(jac.diagonal, 2));
                diagonal(flag, :) = jac.diagonal(N(flag, 1), :);
                
                notFlag = ~flag;
                flag2 = [false(nf, 1); notFlag];
                diagonal(flag2, :) = jac.diagonal(N(notFlag, 2), :);
                jac.diagonal = diagonal;
            else
                diagonal = bsxfun(@times, jac.diagonal(N, :), [flag; ~flag]);
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
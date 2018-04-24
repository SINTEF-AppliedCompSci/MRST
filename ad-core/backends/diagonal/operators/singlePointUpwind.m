function v = singlePointUpwind(flag, N, v)
    % Single-point upwind for the NewAD library
    vD = double(v);
    cells = N(:, 2);
    cells(flag) = N(flag, 1);
    
    value = vD(cells, :);
    if isa(v, 'NewAD')
        M = [];
        DS = [];
        v.val = value;
        for jacNo = 1:numel(v.jac)
            [v.jac{jacNo}, M, DS] = upwindJac(v.jac{jacNo}, flag, N, M, DS);
        end
    else
        v = value;
    end
end

function [jac, M, DS] = upwindJac(jac, flag, N, M, DS)
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
            diagonal = bsxfun(@times, jac.diagonal(N, :), [flag; ~flag]);
            if isempty(jac.subset)
                map = N;
            else
                map = jac.subset(N);
            end
            if isempty(DS)
                jac = DiagonalSubset(diagonal, jac.dim, map);
                DS = jac;
            else
                DS.diagonal = diagonal;
                DS.map = map;
                DS.dim = jac.dim;
%                 DS.nvars = jac.nvars;
                jac = DS;
            end
        end
    end

end
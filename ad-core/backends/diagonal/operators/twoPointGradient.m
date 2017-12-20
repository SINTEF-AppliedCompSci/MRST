function v = twoPointGradient(N, v, M)
% Discrete gradient for the NewAD library
    if isa(v, 'NewAD')
        v.val = gradVal(v.val, N);
        v.jac = cellfun(@(x) gradJac(x, N, M), v.jac, 'UniformOutput', false);
    else
        v = gradVal(v, N);
    end
end

function jac = gradJac(jac, N, M)
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
        if 1
            diagonal = M*jac.diagonal;
        else
            diagonal = jac.diagonal(N, :);
            nf = size(N, 1);
            diagonal(1:nf, :) = -diagonal(1:nf, :);
        end
        if isempty(jac.subset)
            jac = DiagonalSubset(diagonal, jac.dim, N);
        else
            jac = DiagonalSubset(diagonal, jac.dim, jac.subset(N));
        end
    end
end

function v = gradVal(val, N)
    if size(val, 2) > 1
        v = val(N(:, 2), :) - val(N(:, 1), :);
    else
        v = diff(val(N), 1, 2);
    end
end
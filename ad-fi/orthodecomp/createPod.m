function [basis X_bar] = createPod(values, frac_e, n_e)
    [X X_bar] = snapshots(values);
    [V,S,W] = svd(X); % bruk econ?
    eig = diag(S).^2;
    tot_e = sum(eig);
    n = 0;
    energy = 0;
    while true
        n = n + 1;
        energy = energy + eig(n);
%         energy = eig(n);
        if energy/tot_e > frac_e || n == numel(eig) || n == n_e
            break
        end
    end

    fprintf('%d of %d eigenvalues included\n', n, numel(eig));
    basis = V(:, 1:n);
end

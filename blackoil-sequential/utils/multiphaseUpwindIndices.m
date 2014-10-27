function [up, theta, r] = multiphaseUpwindIndices(G, vT, T, K, upstr)
    nPh = numel(G);
    nF  = numel(T);
    G = flattenAndMakeDouble(G);
    K = flattenAndMakeDouble(K);
    nC  = size(K, 1);
    
    [G, sortInd] = sort(G, 2);
    theta = repmat(vT, 1, nPh);
    for i = 1:nPh
        for j = 1:nPh
            if j == i
                continue
            end
            flag = (i < j);
            kj = zeros(nC, 1);
            for k = 1:nPh
                % Upstream weight mobility
                kloc = upstr(flag, K(:, k));
                % Lookup sort index for this phase (G is sorted by weight,
                % K is not)
                ind = sortInd(:, j) == k;
  
                kj(ind) = kloc(ind);
            end
            theta(:, i) = theta(:, i) + T.*(G(:, i) - G(:, j)).*kj;
        end
    end
    [ix, r] = max(theta > 0, [], 2);
    r = r - 1;
    r(all(theta < 0, 2)) = nPh;
    
    up = false(nF, nPh);
    for i = 1:nPh
        up(:, i) = sortInd(:, i) > r;
    end
end

function d = flattenAndMakeDouble(d)
    d = cellfun(@double, d, 'UniformOutput', false);
    d = [d{:}];
end

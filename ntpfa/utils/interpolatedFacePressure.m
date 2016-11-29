function p_f = interpolatedFacePressure(G, p, useDist)
    if nargin == 2
        useDist = true;
    end
    p_set = [0; p];
    n = G.faces.neighbors + 1;
    if useDist
        c = [nan(1, G.griddim); G.cells.centroids];
        
        d_f = zeros(G.faces.num, 1);
        for i = 1:2
            d_f(:, i) = 1./sqrt(sum((c(n(:, i), :) - G.faces.centroids).^2, 2));
        end
        d_f(isnan(d_f)) = 0;
        p_f = bsxfun(@rdivide, sum(p_set(n).*d_f, 2), sum(d_f, 2));
    else
        p_f = sum(p_set(n), 2)./sum(G.faces.neighbors > 0, 2);
    end
end
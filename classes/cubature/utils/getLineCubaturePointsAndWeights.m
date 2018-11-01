function [x, w, n] = getLineCubaturePointsAndWeights(k)
% Get line cubature points and weights for a cubature of precision k

    % Both 0 and 1 uses midpoint rule
    k  = k + (k == 0);
    % Weights calculated from legendre polynomials
    l  = legendrePolynomials(k); 
    dl = cellfun(@(l) dx(l,1), l, 'unif', false);
    dl = dl{end};

    if k <= 1
        x = 0;

    elseif k == 2
        x = [-sqrt(1/3); sqrt(1/3)];

    elseif k == 3
        x = [-sqrt(3/5); 0; sqrt(3/5)];

    elseif k == 4
        x = [-sqrt(3/7 + 2/7*sqrt(6/5)); -sqrt(3/7 - 2/7*sqrt(6/5)); ...
              sqrt(3/7 - 2/7*sqrt(6/5));  sqrt(3/7 + 2/7*sqrt(6/5))];

    elseif k == 5
        a = 5; b = 2*sqrt(10/7);
        x = [-1/3*sqrt(a + b); -1/3*sqrt(a - b); 0; ...
             +1/3*sqrt(a - b); +1/3*sqrt(a + b)];

    elseif k <= 7
        k  = 7;
        l  = legendrePolynomials(k);
        dl = cellfun(@(l) dx(l,1), l, 'unif', false);
        dl = dl{end};
        x = [-0.949107912342758524526189684048;
             -0.741531185599394439863864773281;
             -0.405845151377397166906606412077;
              0;
              0.405845151377397166906606412077;
              0.741531185599394439863864773281;
              0.949107912342758524526189684048];

    end
    
    % Calculate weights
    w = 2./((1 - x.^2).*dl(x).^2);
    w = w/2;
    % Transform coordinates to [0,1]
    x = (x + 1)/2;
    % Assign y coordinates
    x = [x, 1-x];

    n = numel(w);

end
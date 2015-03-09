function T = expandPVCOintoPVTO(tbl, ntab)
% Utility for expanding a PVCO table into a PVTO table
    T = cell([ntab, 1]);
    
    % Relative increments used to evaluate undersaturated curves
    dp = [1.01, 1.05, 1.1, 1.25];
    ndp = numel(dp);
    
    for i = 1:ntab
        tloc = tbl{i};
        npts = size(tloc, 1);
        
        nd = npts*(1 + ndp);
        t = struct('key' , [], ...
                   'pos' , 1, ...
                   'data', zeros(nd, 3));
        t.key = tloc(:, 1);
        
        % We have data for saturated regions, need to evaluate the
        % undersaturated curves
        satix = false(nd, 1);
        satix(1:(ndp+1):(nd - ndp)) = true;
        
        sat_tbl = tloc(:, 2:4);
        
        t.data(satix, :) = sat_tbl;
        t.data(~satix, :) = makeUsatTable(sat_tbl, tloc(:, 5:6), dp);
        t.pos = (1:(1+ndp):(nd + 1))';
        
        T{i} = t;
    end
end
 
function tbl = makeUsatTable(sat_tbl, cv, dp)
    % Compute b/mu for undersaturated curves
    p0 = repmat(sat_tbl(:, 1), 1, numel(dp));
    p = bsxfun(@times, p0, dp);
    
    b  = constantPressureDependence(sat_tbl(:, 2), p, p0, cv(:, 1));
    mu = constantPressureDependence(sat_tbl(:, 3), p, p0, cv(:, 2));
    
    fix = @(x) reshape(x', [], 1);
    tbl = [fix(p), fix(b), fix(mu)];
end

function v = constantPressureDependence(v0, p, p0, c)
    % v(p) = v(p_0) e ^ -(p - p_0) C
    pdiff = p0 - p;
    e = exp(bsxfun(@times, pdiff, c));
    v = bsxfun(@times, v0, e);
end
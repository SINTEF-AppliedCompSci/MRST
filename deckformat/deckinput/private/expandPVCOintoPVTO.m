function T = expandPVCOintoPVTO(tbl, ntab)
% Utility for expanding a PVCO table into a PVTO table
    T = cell([ntab, 1]);
    
    for i = 1:ntab
        tloc = tbl{i};
        npts = size(tloc, 1);
        
        t = struct('key' , [], ...
                   'pos' , 1, ...
                   'data', zeros(npts*2, 3));
        t.key = tloc(:, 1);
        
        % We have data for saturated regions, need to get interpolated
        % values for undersaturated retiongs
        satix = 1:2:(2*npts-1);
        usatix = 2:2:(2*npts);
        
        sat_tbl = tloc(:, 2:4);
        
        t.data(satix, :) = sat_tbl;
        t.data(usatix, :) = makeUsatTable(sat_tbl, tloc(:, 5:6));
        t.pos = (1:2:(2*npts + 1))';
        
        T{i} = t;
    end
end
 
function tbl = makeUsatTable(sat_tbl, cv)
    p0 = sat_tbl(:, 1);
    p = 1.25*p0;
    
    % Interpolate using constants
    b  = sat_tbl(:, 2).*exp(-(p - p0).*cv(:, 1));
    mu = sat_tbl(:, 3).*exp(-(p - p0).*cv(:, 2));
    
    tbl = [p, b, mu];
end
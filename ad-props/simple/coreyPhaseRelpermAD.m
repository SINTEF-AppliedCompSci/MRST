function fn = coreyPhaseRelpermAD(n, sr, kwm, sr_tot)
% Make function handle for AD style fluids

    % n - exponent
    % sr - residual saturation
    % kwm - endpoint relperm
    % sr_tot - sum of sr for all phases present
    
    if nargin < 1
        n = 1;
    end
    if nargin < 2
        sr = 0;
    end
    if nargin < 3
        kwm = 1;
    end
    if nargin < 4
        sr_tot = sr;
    end
%     den = 1 - sr_tot;
%     fn = @(s) kwm*max(min(((s - sr)./den), 1), 0).^n;
    fn = @(s) coreyRelperm(s, n, sr, kwm, sr_tot);
end

function kr = coreyRelperm(s, n, sr, kwm, sr_tot)
    den = 1 - sr_tot;
    sat = ((s - sr)./den);
    if isa(sat, 'ADI')
        sat.val = max(min(sat.val, 1), 0);
    else
        sat = max(min(sat, 1), 0);
    end
    
    kr = kwm*sat.^n;
end

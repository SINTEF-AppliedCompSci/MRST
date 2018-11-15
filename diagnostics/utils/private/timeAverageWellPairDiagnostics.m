function WP = timeAverageWellPairDiagnostics(d, timeSteps)
assert(numel(timeSteps)>1);
if nargin < 2
    timeSteps = (1:numel(d.Data.diagnostics))';
end

t  = d.Data.time.cur(timeSteps);
dt = diff(t(:));
w  = dt/sum(dt);
diagn = d.Data.diagnostics(timeSteps(2:end));
WP    = diagn(1).WP;
for k = 1:(numel(timeSteps)-1)
    z = k>1;
    for ni = 1:numel(WP.inj)
        WP.inj(ni).alloc  = z*WP.inj(ni).alloc  + w(k)*diagn(k).WP.inj(ni).alloc;
        WP.inj(ni).ralloc = z*WP.inj(ni).ralloc + w(k)*diagn(k).WP.inj(ni).ralloc;
    end
    
    for np = 1:numel(WP.prod)
        WP.prod(np).alloc  = z*WP.prod(np).alloc  + w(k)*diagn(k).WP.prod(np).alloc;
        WP.prod(np).ralloc = z*WP.prod(np).ralloc + w(k)*diagn(k).WP.prod(np).ralloc;
    end
    
    WP.vols = z*WP.vols + w(k)*diagn(k).WP.vols;
end

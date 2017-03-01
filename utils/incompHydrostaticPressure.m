function [p, fn] = incompHydrostaticPressure(G, contacts, densities, varargin)
    opt = struct('topPressure', 0, ...
                 'nsamples', 1000);
    opt = merge_options(opt, varargin{:});
    ns = opt.nsamples;
    pdist = zeros(ns, 1);
    zi = zeros(ns, 1);
    zmax = max(G.cells.centroids(:, 3));
    zmin = min(G.cells.centroids(:, 3));
    h = (zmax - zmin)./(ns-1);
    pdist(1) = opt.topPressure;
    zi(1) = zmin;
    for i = 2:ns
        zi(i) = h*i + zmin;
        ix = find(contacts < zi(i), 1, 'last');
        if isempty(ix)
            ix = 0;
        end
        rho = densities(ix+1);

        pdist(i) = pdist(i-1) + rho*h*norm(gravity);
    end
   
    z = G.cells.centroids(:, 3);
    
    fn = @(z) interp1(zi, pdist, z, 'linear', 'extrap');
    p = fn(z);
end
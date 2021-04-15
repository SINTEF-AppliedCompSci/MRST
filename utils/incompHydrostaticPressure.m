function [p, fn] = incompHydrostaticPressure(G, contacts, densities, varargin)
%Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

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

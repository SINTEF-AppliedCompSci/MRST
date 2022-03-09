function tab = postprocessTable(tab, varargin)
%Postprocess table created using generateCoolPropTable

    opt = struct('fixPatch' , false, ...
                 'rangeX'   , []   , ...
                 'rangeY'   , []   , ...
                 'primaryAx', 'x'  , ...
                 'fixZeros' , true , ...
                 'smoothen' , true );
    opt = merge_options(opt, varargin{:});
    
    if opt.fixPatch
        tab = fixPatch(tab, opt);
    end
    if opt.fixZeros
        tab = fixZeros(tab, opt);
    end
    if opt.smoothen
        tab = smoothen(tab, opt);
    end
    
end

%-------------------------------------------------------------------------%
function tab = fixPatch(tab, opt)
    v = opt.rangeY(1);
    ixX = tab.x >= opt.rangeX(1) & tab.x <= opt.rangeX(2);
    [~, ixY] = min(abs(tab.y - v));
    while ixY <= numel(tab.y)
        z = tab.data(~ixX, ixY);
        tab.data(ixX, ixY) = interp1(tab.x(~ixX), z, tab.x(ixX));
        ixY = ixY + 1;
    end
end

%-------------------------------------------------------------------------%
function tab = fixZeros(tab, opt)
    n = numel(tab.x);
    data = tab.data;
    fix  = data == 0;
    if ~any(any(fix)), return; end
    [x,y] = ndgrid(tab.x, tab.y);
    fn = scatteredInterpolant(x(~fix), y(~fix), data(~fix), 'linear', 'nearest');
    data = fn(x, y);
    tab.data = reshape(data, n, n);
end

%-------------------------------------------------------------------------%
function tab = smoothen(tab, opt)
    n = numel(tab.x);
    G = cartGrid([n,n]);
    N = getNeighbourship(G);
    N = [N; fliplr(N)];
    nc = accumarray(N(:,1), 1);
    h = 1./n;
    S =  spdiags(nc, 0, n^2, n^2) - sparse(N(:,1), N(:,2), 1, n^2, n^2);
    data = tab.data(:);
    
    d = data./max(data);
    dd = (S*d)./h^2;
    fix = abs(dd) > 1e4;
    
    if ~any(fix), return; end
    data0 = data;
    nc = accumarray(N(:,1), ~fix(N(:,2)));
    M = sparse(N(:,1), N(:,2), ~fix(N(:,2)), n^2, n^2);
    data = M*data;
    data = data./nc;
    ix = nc == 0;
    data(ix) = data0(ix);
    tab.data = reshape(data, n, n);
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
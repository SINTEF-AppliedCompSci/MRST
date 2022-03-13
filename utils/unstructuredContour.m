function [c, h] = unstructuredContour(G, val, nl, varargin)
%Undocumented Utility Function

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

    opt = struct('useNodeMax', false   , ...
                 'fill'      , false   , ...
                 'n'         , [20, 20], ...
                 'extrapolationMethod', 'none');
    [opt, extra] = merge_options(opt, varargin{:});

    if nargin == 2
        nl = 10;
    end
    fun = scatteredInterpolant(G.cells.centroids(:,1), G.cells.centroids(:,2), val, 'linear', opt.extrapolationMethod);
    if opt.useNodeMax
        xMax = max(G.nodes.coords);
        xMin = min(G.nodes.coords); 
    else
        xMax = max(G.cells.centroids);
        xMin = min(G.cells.centroids);
    end
    [x, y] = ndgrid(linspace(xMin(1), xMax(1), opt.n(1)), linspace(xMin(2), xMax(2), opt.n(2)));
    val = fun(x, y);
    
    if opt.fill
        [c, h] = contourf(x, y, val, nl, extra{:});
    else
        [c, h] = contour(x, y, val, nl, extra{:});
    end
    
end

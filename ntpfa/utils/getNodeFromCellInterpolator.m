function [P, nodeNo, cellNo, w] = getNodeFromCellInterpolator(G)
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

    % Get nodes
    [nodeNo, pos] = gridCellNodes(G, 1:G.cells.num, 'unique', true);
    cellNo = rldecode((1:G.cells.num)', diff(pos));
    
    % One over distance weighting
    w = 1./sqrt(sum((G.cells.centroids(cellNo, :) - G.nodes.coords(nodeNo, :)).^2, 2));
    
    % Divide by sum of weights for partition of unity
    sumw = accumarray(nodeNo, w);
    w = w./sumw(nodeNo);
    
    % n_n by n_c matrix constructing values on nodes by cells
    P = sparse(nodeNo, cellNo, w, G.nodes.num, G.cells.num);
end

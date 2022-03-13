function [F, mask_fn] = continuousCelldataFunction(Gt, celldata, quickclip)
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

    % computing mask function
    bnd = Gt.nodes.coords(extractGridBoundaryNodes(Gt),:);
    mask_fn = @(x, y) inpolygon(x, y, bnd(:,1), bnd(:,2));

    % computing interpolating function
    if ~quickclip
        % we will not use boundary notes as NaN values for 'quick' clipping
        bnd = zeros(0,2); 
    end
    
    x = [Gt.cells.centroids(:, 1); bnd(1:end-1,1)];
    y = [Gt.cells.centroids(:, 2); bnd(1:end-1,2)];
    z = [celldata; ones(size(bnd,1)-1,1) * NaN];
    
    F = TriScatteredInterp(x, y, z);
end

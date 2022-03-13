function colorizeRegion(h, Gt, field, color, varargin)
% Colorize area described by positive value for field

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

    samples = 200;
    figure(h);
    
    field   = double( field > 0);
    gx      = Gt.cells.centroids(:,1);
    gy      = Gt.cells.centroids(:,2);
    % following lines creates a band of NaN values around grid boundary, to
    % prevent extrapolation outside in the case of concave shapes.
    n       = extractGridBoundaryNodes(Gt);
    n       = n(1:end-1);
    corners = Gt.nodes.coords(n,:);
    % Create interpolant
    F       = TriScatteredInterp([gx(:); corners(:,1)], ...
                                 [gy(:); corners(:,2)], ...
                                 [field(:); NaN*ones(size(n))] );
    % Sampling grid
    X       = linspace(min(gx), max(gx), samples);
    Y       = linspace(min(gy), max(gy), samples);
    [x, y]  = meshgrid(linspace(min(gx), max(gx), samples), ...
                       linspace(min(gy), max(gy), samples));
    % Computing region contour
    cont    = contourc(X, Y, F(x,y), [0.1, 0.1]);
    
    i = 1;
    while i < size(cont, 2)
        num_points = cont(2,i);
        xp = cont(1, (i+1):(i+num_points));
        yp = cont(2, (i+1):(i+num_points));
        patch(xp, yp, color, varargin{:});
        i = i + num_points + 1;
    end
end

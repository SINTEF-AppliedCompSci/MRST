function [hs, htop] = plotCurveTube(curve, radius, varargin)
%Plot a curve as a tube with variable width
%
% SYNOPSIS:
%   plotCurveTube(points, radius)
%   plotCurveTube(points, radius, 'data', data)
%
% DESCRIPTION:
%   Plot a curve given by consecutive points as a variable radius tube,
%   with optional colorization based on some underlying dataset. Typical
%   usage is for plotting of well trajectories.
%
% REQUIRED PARAMETERS:
%   curve   - N by 3 array where each row corresponds to the 3D coordinates
%             of points along the curve.
%
%   radius  - Either one value per point in "curve", or a single point
%             to be used for all points. The radius determines the radius
%             of the tube plotted.
%
% OPTIONAL PARAMETERS:
%   data    - Dataset with one value per curve point given as a column
%             vector. Will be used to colorize the tube based on the
%             current figure colormap.
%
%   interpMethod - Type of interpolation used when refining the curve. Uses
%                  standard Matlab interp1 options. Default: Linear.
%
%   interpDataMethod - Type of interpolation used to interpolate data
%                      between points. Default: Linear.
%
%   EdgeColor - Edge color for final patch object.
%
%   Refine    - Refinement factor for curve. Default to 1, which is no
%               refinement. Larger values will give a smoother curve when
%               using e.g. spline interpolation.
%
% RETURNS:
%   Nothing.
%
% EXAMPLE:
%    pts = [0, 0, 5; 0 0 -1; .25 .25 -5; .25, .25, -10];
%    plotCurveTube(pts, 0.1, 'data', pts(:, 3), 'interpMethod', 'spline',...
%                 'refine', 10)
%
% SEE ALSO:
%   

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
    opt = struct('data',    [], ...
                 'interpMethod', 'linear', ...
                 'interpDataMethod', 'linear', ...
                 'EdgeColor', {'none'}, ...
                 'color',   {[1 0 0]}, ...
                 'refine',  1);
    opt = merge_options(opt, varargin{:});
    
    [x, y, z] = cylinder(1);
    C = [z(1, :)', y(1, :)', x(1, :)'];
    
    nel = size(curve, 1);

    
    if size(radius, 1) == 1 && nel ~= 1
        radius = repmat(radius, nel, 1);
    end
    
    data = opt.data;
    if opt.refine ~= 1 && nel > 1
        d = cumsum([0; sqrt(sum(diff(curve).^2, 2))]);
        len = d(end);
        d_new = (0:len/(nel*opt.refine):len)';
        % Include old data points
        d_new = unique([d_new; d]);
        d_new = sort(d_new);

        curve  = interp1(d, curve,  d_new, opt.interpMethod);
        radius = interp1(d, radius, d_new, opt.interpMethod);
        if ~isempty(data)
            data = interp1(d, data, d_new, opt.interpDataMethod);
        end
        nel = size(curve, 1);
    end
    
    if size(curve, 2) == 2
        curve = [curve, zeros(nel, 1)];
        radius = [radius, mean(radius, 2)];
    end
    
    hasData = ~isempty(data);
    
    if nel == 1
        % We were only given a single data point - extend it to two almost
        % identical data points
        vectors = [0 0 1; 0 0 1];
        curve = [curve; curve];
        curve(:, 3) = curve(:, 3) +  sqrt(eps)*[-1; 1];
        if hasData
            data = [data; data];
        end
        radius = [radius; radius];
        nel = 2;
    else
        % Tangential vectors along the curve using first order
        % approximation
        vectors = diff(curve, 1, 1);
        % First element repeated to ensure that we have a vector for each
        % data point
        vectors = [vectors(1, :); vectors];
        vectors = normalize(vectors);
    end
    [x, y, z] = deal(zeros(nel, size(C, 1)));
    
    for i = 1:nel
        % Loop through all the segments and rotate the reference circles,
        % adding them to the list of cylinders as we go
        v = vectors(i, :);
        % Rotate
        pts = rotate(C, v);
        pts = bsxfun(@mtimes, pts, radius(i, :));
        % Translate
        pts = bsxfun(@plus, pts, curve(i, :));

        sub = i;
        x(sub, :) = pts(:, 1)';
        y(sub, :) = pts(:, 2)';
        z(sub, :) = pts(:, 3)';
    end
    
    if hasData
        d = repmat(data, 1, size(z, 2));
        hs = surf(x, y, z, d, 'EdgeColor', opt.EdgeColor);
    else
        hs = surf(x, y, z, 'EdgeColor', opt.EdgeColor, 'FaceColor', opt.color, 'CData', []);
    end
    
    % Plot circles at the endpoints
    for i = 1 %[1, size(x, 1)];
        if hasData
            htop = patch(x(i, :), y(i, :), z(i, :), d(i, :));
        else
            htop = patch(x(i, :), y(i, :), z(i, :), nan, 'FaceColor', opt.color);
        end
    end
end
        
function pts = rotate(pts, newvec)
    angle = acos(newvec);
    M = rotateX(angle(1))*rotateY(angle(2));
    for i = 1:size(pts, 1)
        pts(i, :) = pts(i, :)*M;
    end
end

function M = rotateX(o)
    M = [1,      0,        0; ...
         0, cos(o),  -sin(o); ...
         0, sin(o),   cos(o);];
end

function M = rotateY(o)
    M = [cos(o),  0, sin(o); ...
         0,     1,     0; ...
         -sin(o), 0, cos(o);];
end

function M = rotateZ(o)
    M = [cos(o),  -sin(o), 0; ...
         sin(o), cos(o), 0; ...
         0, 0, 1;];
end

function vectors = normalize(vectors)
    vectors = bsxfun(@rdivide, vectors, sqrt(sum(vectors.^2, 2)));
end

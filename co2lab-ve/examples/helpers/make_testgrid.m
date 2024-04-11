function G = make_testgrid(res, extent, bend, slope, depth, bumps, varargin)
% Make a simple test grid, sloping along the x-direction and with a slight y
% curvature (funnel-shaped).  Bumps may also be added, on the format [xpos,
% width, height], where 'xpos' is relative along the length of the grid (from 0
% to 1), 'width' is the width of the trap (relative to the grid size in the
% y-dimension), and 'heigh' is the amplitude of the trap (in meters)

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

    opt.ripples_amplitude_x = 0;
    opt.ripples_number_x = 100;
    opt = merge_options(opt, varargin{:});
    
    %% Creating base sloping surface with slight y curvature
        
    [X, Y, Z] = deal(extent(1), extent(2), extent(3));
    [resx, resy, resz] = deal(res(1), res(2), res(3));
    
    
    x = linspace(0, X, resx + 1);
    y = linspace(0, Y, resy + 1);
    
    base_z = depth - (x * slope)' .* ones(1, numel(y));
    ybend = bend * (y/Y) .* (1 - y/Y);
    base_z = bsxfun(@minus, base_z, ybend);
    
    %% Adding traps
    
    [xx, yy] = meshgrid(x, y);
    xx = xx';
    yy = yy';
    
    z = base_z;    
    
    for b = bumps'
    
        xpos = b(1) * X;
        ypos = 0.5 * Y;
        width = b(2) * Y;
        height = b(3);

        
        dist_x = xx - xpos;
        if res(2) == 1
            % cross-sectional problem
            dist_y = 0;
        else
            dist_y = yy - ypos;
        end
        
        d2 = dist_x.^2 + dist_y.^2;
        
        trap = height * exp(-d2/width^2);
        
        z = z - trap;
    end
    
    
    %% Constructing the mrst grid
    G = cartGrid(res, extent);
    G.nodes.coords(:,3) = G.nodes.coords(:,3) + repmat(z(:), resz+1, 1);
    
    %% Add ripples if requested
    if opt.ripples_amplitude_x > 0
        G = add_topsurface_ripples(G, opt.ripples_amplitude_x, opt.ripples_number_x);
    end

    %% Compute final geometry and return
    G = computeGeometry(G);

end

% ----------------------------------------------------------------------------
function G = add_topsurface_ripples(G, amplitude, number)
    
    N = G.cartDims(1) + 1; % number of corners along x-axis in top surface
    
    %% computing top surface ripples
    ripples = amplitude * sin(linspace(0, 1, N) * 2 * pi * number);
    
    %% adjusting grid nodes
    dz = repmat((linspace(1, 0, G.cartDims(3)+1).^2)' * ripples, 1, 2)';
    
    G.nodes.coords(:,3) = G.nodes.coords(:,3) + dz(:);
end

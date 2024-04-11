function G = make_testgrid(res, extent, bend, slope, depth, bumps)

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
    
    G = computeGeometry(G);

end

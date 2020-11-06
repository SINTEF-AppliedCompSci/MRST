function G = createEquilateralTriangleGrid(Nx)

    dx = 1/Nx(1); 
    x_size = 1; 
    y_size = 1; 
    x_min = 0; y_min = 0; 

    dy = sqrt (3) / 2 * dx; 
    % adjust so that we get a round nice number of points in each direction
    nx = max (1, round (x_size / dx)); 
    dx = x_size / nx; 
    ny = max (1, round (y_size / dy)); 
    dy = y_size / ny; 

    % write vertices at these points; we have one more point than there are
    % sides, which is we we start at zero instead of one
    vertices = []; 
    for i = 0:ny
        % every other row is adjusted a half column to the right and the last
        % point disappears.
        if ~mod(i, 2)
            ofs = 0; 
            count = nx; 
        else
            ofs = dx / 2; 
            
            % if we only adjust the row, we'll loose two point along each edge
            % and this will cause a long edge at the boundary. as we prefer
            % shorter ones, we'll setup two points at the edges as well
            vertices = [vertices; x_min, y_min + i*dy; x_min + nx*dx   y_min + i*dy]; 
            
            % since we are adjusting with the offset; don't write the last end
            % point (thus minus one on the upper bound). we handle the end
            % point of the line above.
            count = nx - 1; 
        end; 
        
        % disperse the points along the line.
        ind = (0 : count)'; 
        vertices = [ vertices; x_min + ind * dx + ofs, y_min + i * ones(count + 1, 1) * dy]; 
        
    end; 
    p = vertices; 
    t = delaunayn(p); 
    g = triangleGrid(p, t); 
    G = computeGeometry(g); 
end
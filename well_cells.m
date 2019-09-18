function cells = well_cells(G, type)

    mrstModule add wellpaths
    
% seems like first node is off for eni grid?

    xyzmin = min(G.nodes.coords(2:end,:));
    xyzmax = max(G.nodes.coords(2:end,:));
    xmin = xyzmin(1);
    ymin = xyzmin(2);
    zmin = xyzmin(3);
    xmax = xyzmax(1);
    ymax = xyzmax(2);
    zmax = xyzmax(3);
    
    volume_hmin = min(G.cells.volumes.^(1/3));
    area_hmin = min(sqrt(G.faces.areas));
    hmin = min(volume_hmin, area_hmin);
    
    % offset differently in x and y
    offset = 1e-1*hmin;

    switch type
        case '1' % at xy min, in z direction
            xyz0 = xyzmin + offset;
            xyz1 = [xmin+2*offset, ymin+offset, zmax-offset];
            
        case '2' % in y direction, located at (xmax, zmin)
            xyz0 = [xmax-2*offset, ymin+offset, zmin+offset];
            xyz1 = [xmax-2*offset, ymax-offset, zmin+offset];
            
        case '3' % center (minus offset to avoid grid alignment)
            xyz0 = [0.5*(xmin+xmax)-offset, 0.5*(ymin+ymax)-offset, zmin+offset];
            xyz1 = [0.5*(xmin+xmax)-offset, 0.5*(ymin+ymax)-offset, zmax-offset];
            
        case '4' % at x max, ymin (matches 1)
            xyz0 = [xmax-2*offset, ymin+offset, zmin+offset];
            xyz1 = [xmax-2*offset, ymin+offset, zmax-offset];
            
      case '5' % diagonally across 1
        xyz0 = [xmax-offset, ymax-offset, zmin+offset];
        xyz1 = [xmax-offset, ymax-offset, zmax-offset];
            
        otherwise
            disp('unknown case');
        keyboard;
    end
        
    wellpath = makeSingleWellpath([xyz0; xyz1]);
    cells = findWellPathCells(G, wellpath);
    
    % % Without wellpaths:
    % points_per_cell = 2;
    % N = points_per_cell*round((zmax-zmin)/hmin);
    % p = xyz0 + linspace(0,1,N).'.*(xyz1-xyz0);
    
    % cells = zeros(N, 1);
    % for i = 1:N
    %     cells(i) = findEnclosingCell(G, p(i,:));
    % end
    % cells = unique(cells)';

    % % debug
    % figure
    % plotGrid(G,'facealpha',0.2);
    % hold on
    % plot3(p(:,1),p(:,2),p(:,3),'.');
    % keyboard
    
end

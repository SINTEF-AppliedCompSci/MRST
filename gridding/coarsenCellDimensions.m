function G = coarsenCellDimensions(G, varargin)

    opt = struct('scaleDx'    , false, ...
                 'useBBcoords', false);
    opt = merge_options(opt, varargin{:});

    % Get minimum and maximum cell coordinates
    [~, ix]  = sort(G.partition);
    [~ , nc] = rlencode(G.partition(ix),1);
    [G.cells.xMin, ~           ] = getMinMax(G.parent.cells.xMin(ix,:), nc);
    [~           , G.cells.xMax] = getMinMax(G.parent.cells.xMax(ix,:), nc);
    % Calculate bounding box dimensions
    G.cells.dx = G.cells.xMax - G.cells.xMin;
    
    % Get minimum and maximum face coordinates
    nf = diff(G.faces.connPos);
    [G.faces.xMin, ~           ] = getMinMax(G.parent.faces.xMin(G.faces.fconn,:), nf);
    [~           , G.faces.xMax] = getMinMax(G.parent.faces.xMax(G.faces.fconn,:), nf);
    % Calculate bounding box dimensions
    G.faces.dx = G.faces.xMax - G.faces.xMin;
    
    n = G.faces.normals./sqrt(sum(G.faces.normals.^2,2));
    if G.griddim == 2
        v = [n(:,2), -n(:,1)];
        x = [G.faces.centroids - G.faces.areas/2.*v;
             G.faces.centroids + G.faces.areas/2.*v];
        ix = [1:G.faces.num; (1:G.faces.num) + G.faces.num];
        x = x(ix(:),:);
    else
        if ~any(strcmpi(G.type, 'generateHybridCoarseGrid'))
            nodes = G.parent.faces.nodes(mcolon(G.parent.faces.nodePos(G.faces.fconn),...
                                                G.parent.faces.nodePos(G.faces.fconn+1)-1));

            nfn = diff(G.parent.faces.nodePos);
            ix  = rldecode((1:G.faces.num)', diff(G.faces.connPos), 1);
            nfn = accumarray(ix, nfn(G.faces.fconn));

            % Construct auxilliary coordinate system
            vec1 = [n(:,2) - n(:,3), -(n(:,1) + n(:,3)), n(:,1) + n(:,2)];
            vec1 = vec1./sqrt(sum(vec1.^2,2));
            n    = G.faces.normals./sqrt(sum(G.faces.normals.^2,2));
            vec2 = cross(vec1, n, 2);

            xn = G.parent.nodes.coords(nodes,:);
            xn = xn - rldecode(G.faces.centroids, nfn, 1);
            xn = [sum(xn.*rldecode(vec1, nfn, 1),2), sum(xn.*rldecode(vec2, nfn, 1),2)];

            % Construct new coordinate system aligned with vector from center
            % to maximum coordinate (~approximate minimum bounding rectangle)        
            [xMin, xMax] = getMinMax(xn, nfn);
            xc = (xMax - xMin)/2;
            vec1 = xc(:,1).*vec1 + xc(:,2).*vec2;
            vec1 = vec1./sqrt(sum(vec1.^2,2));
            vec2 = cross(vec1, n, 2);

            % Save coordinate system
            G.faces.coordSys = {vec1, vec2};

            % Calculate face node coordinates
            xn = G.parent.nodes.coords(nodes,:);
            xn = xn - rldecode(G.faces.centroids, nfn, 1);
            xn = [sum(xn.*rldecode(vec1, nfn, 1),2), sum(xn.*rldecode(vec2, nfn, 1),2)];
            [xMin, xMax] = getMinMax(xn, nfn);
            dx = xMax - xMin;
            if opt.scaleDx
                % Optional scaling of dx so that prod(dx,2) = G.faces.areas
                dx = dx.*sqrt(G.faces.areas./prod(dx,2));
            end
            G.faces.dx = dx;

            if opt.useBBcoords
                x = [-dx(:,1)/2.*vec1 - dx(:,2)/2.*vec2; 
                      dx(:,1)/2.*vec1 - dx(:,2)/2.*vec2;
                      dx(:,1)/2.*vec1 + dx(:,2)/2.*vec2;
                     -dx(:,1)/2.*vec1 + dx(:,2)/2.*vec2];
            else
                x = [xMin(:,1), xMin(:,2);
                     xMax(:,1), xMin(:,2);
                     xMax(:,1), xMax(:,2);
                     xMin(:,1), xMax(:,2)];
                x = x(:,1).*repmat(vec1,4,1) + x(:,2).*repmat(vec2,4,1);
            end

            % Translate to face centroids
            x  = x + repmat(G.faces.centroids, 4, 1);
            ix = repmat(1:G.faces.num, 4, 1) + (0:3)'*G.faces.num;
            x  = x(ix(:), :);
        end
    end
    
    % Assign face nodes
    
    if ~any(strcmpi(G.type, 'generateHybridCoarseGrid'))
        G.nodes.coords = x;
        G.nodes.num    = size(x,1);
        G.faces.nodes  = (1:size(x,1))';
        stp = 2*(G.griddim-1);
        pos = 1:stp:stp*(G.faces.num+1);
        G.faces.nodePos = pos;
        
        if G.griddim == 3
            
            nodes = [1:size(G.nodes.coords,1), 1:G.faces.num]';
            G.edges.nodes   = nodes;
            G.edges.num     = numel(nodes)/2;
            G.edges.nodePos = (1:2:2*(G.edges.num+1))';
            G.faces.edges   = (1:G.edges.num)';
            G.faces.edgePos = (1:4:4*(G.faces.num+1))';
            
        end
        
    end
    
end

function coordSys = faceCoordSys(G)
    
    % Get nodes of last edge of each face
    nfe   = diff(G.faces.edgePos);
    edges = G.faces.edges(cumsum(nfe));
    nodes = G.edges.nodes(mcolon(G.edges.nodePos(edges), ...
                                 G.edges.nodePos(edges+1)-1));

    % First coordinate axis points along last edge
    vec1 = G.nodes.coords(nodes(2:2:end),:) ...
         - G.nodes.coords(nodes(1:2:end-1),:);
    vec1 = vec1./sqrt(sum(vec1.^2,2));
    
    % Second vector constructed as cross product of vec1 and face normal
    n    = G.faces.normals./G.faces.areas;
    vec2 = cross(vec1, n, 2);
    
    coordSys = {vec1, vec2};

end

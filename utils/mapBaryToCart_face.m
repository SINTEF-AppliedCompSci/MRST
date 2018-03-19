function [b2c, nSimplex, vol] = mapBaryToCart_face(G, cells, type)

    faces = G.cells.faces(mcolon(G.cells.facePos(cells), G.cells.facePos(cells+1)-1));

    switch type
        
        case 'volume'
            
%             faces = G.cells.faces(:,1);
            if G.griddim == 2
                
                nodes = G.faces.nodes(mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces+1)-1));
                
                % One triangle per face of each cell
                nSimplex = diff(G.cells.facePos);
                nSimplex = nSimplex(cells);
                
            else
                edges = G.faces.edges(mcolon(G.faces.edgePos(faces), G.faces.edgePos(faces+1)-1));
                nodes = G.edges.nodes(mcolon(G.edges.nodePos(edges), G.edges.nodePos(edges+1)-1));
                nfe = diff(G.faces.edgePos);
                
                % One triangle per edge of each face of each cell
                nSimplex = accumarray(rldecode((1:G.cells.num)', diff(G.cells.facePos), 1), nfe(G.cells.faces(:,1)));
                nSimplex = nSimplex(cells);
                
            end
            
            
            xn = G.nodes.coords(nodes,:);
            xf = [];
            xc = rldecode(G.cells.centroids(cells,:), nSimplex, 1);
            
            
            nspts = G.griddim + 1; % Number of simplex points (2 for lines, 3 for triangles, 4 for tetrahedra)
            
            vec = 1:nspts:sum(nSimplex)*nspts-nspts+1;
            ixn = mcolon(vec, vec+1);
            ixf = [];
            ixc = vec + nspts - 1;

            if G.griddim == 3
                xf  = rldecode(G.faces.centroids(faces,:), nfe(faces), 1);
                ixf = vec + nspts - 2;
            end

            x([ixn, ixf, ixc],:) = [xn; xf; xc];

            R = zeros(nspts, G.griddim*sum(nSimplex));
            for dNo = 1:G.griddim
                R(:,dNo:(G.griddim):end) = reshape(x(:,dNo), nspts, []);
            end

            b2c = @(xb) map(G, R, xb);

            x0 = x(vec,:);
            v = cell(G.griddim,1);
            for dNo = 1:G.griddim
                v{dNo} = x0 - x(vec + dNo, :);
            end

            if G.griddim == 2
                vol = abs(v{1}(:,1).*v{2}(:,2) - v{1}(:,2).*v{2}(:,1))/2;
            else
                vol = abs(dot(cross(v{1}, v{2}, 2), v{3}, 2))/6;
            end
            
        case 'surface'
            
            % Boundary faces are not included
            [bf, bc] = boundaryFaces(G);
            ncbf     = sum((1:G.cells.num)' == bc',2);
            ncf      = diff(G.cells.facePos);
            ncf      = ncf(cells) - ncbf(cells);
            faces    = faces(~ismember(faces, bf));
            
            
            if G.griddim == 2
            
                nodes = G.faces.nodes(mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces+1)-1));
                xn    = G.nodes.coords(nodes,:);
                
                R = zeros(G.griddim, G.griddim*sum(ncf));
                for dNo = 1:G.griddim
                    R(:, dNo:G.griddim:end) = reshape(xn(:,dNo), G.griddim, []);
                end

                nSimplex = ones(numel(faces),1);
    
            else
                
                edges = G.faces.edges(mcolon(G.faces.edgePos(faces), G.faces.edgePos(faces+1)-1));
                nodes = G.edges.nodes(mcolon(G.edges.nodePos(edges), G.edges.nodePos(edges+1)-1));
                nSimplex = diff(G.faces.edgePos);
                nSimplex = nSimplex(faces);
                
            end
            
            xn = G.nodes.coords(nodes, :);
            xf = [];
            
            nspts = G.griddim;
            vec = 1:nspts:sum(nSimplex)*nspts-nspts+1;
            ixn = mcolon(vec, vec+1);
            ixf = [];

            if G.griddim == 3
                xf = rldecode(G.faces.centroids(faces,:), nSimplex, 1);
                ixf = vec + nspts - 1;
            end

            x([ixn, ixf],:) = [xn; xf];

            R = zeros(nspts, G.griddim*sum(nSimplex));
            for dNo = 1:G.griddim
                R(:, dNo:G.griddim:end) = reshape(x(:,dNo), nspts, []);
            end
            
            b2c = @(xb) map(G, R, xb);
                
            if G.griddim == 2
                
                vol = G.faces.areas(faces);
                
            else
                
                x0 = x(vec,:);
                v = cell(G.griddim,1);
                for dNo = 1:G.griddim-1
                    v{dNo} = x0 - x(vec + dNo, :);
                end

                vol = sqrt(sum(cross(v{1}, v{2}, 2).^2,2))/2;
                
            end

    end    
    
end
    
function x = map(G, R, xb)

     xx = xb*R;
     
     x = zeros(numel(xx)/G.griddim, G.griddim);
     for dNo = 1:G.griddim
        xtmp = xx(:, dNo:G.griddim:end);
        x(:,dNo) = xtmp(:);
     end
          
end
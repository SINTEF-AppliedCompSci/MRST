classdef TriangleCubature < Cubature
    % Triangle cubature class for cubatures on MRST grids
    
    properties
        
        triangulation   % Struct containing triangulation information
        
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function cubature = TriangleCubature(G, prescision)
            % Set up triangle cubature
            cubature = cubature@Cubature(G, prescision, 2);
        end
        
        %-----------------------------------------------------------------%
        function tri = getTriangulation(cubature, G) %#ok
            % Triangulate all cells or faces of G
            
            if G.griddim == 2
                % This is a cell cubature
                % All faces of all cells
                faces   = G.cells.faces(:, 1);
                % All nodes of all faces of all cells
                nodes   = G.faces.nodes(mcolon(G.faces.nodePos(faces)     , ...
                                               G.faces.nodePos(faces+1)-1));
                % Number of triangles per cell
                nTri    = diff(G.cells.facePos);
                parents = rldecode((1:G.cells.num)', nTri, 1);
                % Cell centroids
                xc      = G.cells.centroids(parents, :);
            else
                % This is a face cubature
                % All edges of all faces
                edges   = G.faces.edges;
                % All nodes of all edges of all faces
                nodes   = G.edges.nodes(mcolon(G.edges.nodePos(edges), ...
                                               G.edges.nodePos(edges+1)-1));
                % Number of triangles per face
                nTri    = diff(G.faces.edgePos);
                parents = rldecode((1:G.faces.num)', nTri, 1);
                % Face centroids
                xc      = G.faces.centroids(parents, :);
            end

            % All vertices of triangulation
            xn     = G.nodes.coords(nodes, :);
            coords = [xn; xc];
            % Number of vertices for a triangle
            nPts       = 3; 
            numNodes   = numel(nodes);
            numParents = numel(parents);
            % Index into coordinates for each triangle
            vertPos     = zeros(sum(nTri)*nPts, 1);
            % First two coordinates of each triangle are nodes
            vec         = 1:nPts:sum(nTri)*nPts;
            ii          = mcolon(vec, vec+1);
            vertPos(ii) = 1:numNodes;
            % Last cordinate is the cell centroid
            ii          = nPts:nPts:sum(nTri)*nPts;
            vertPos(ii) = (1:numParents) + numNodes;
            % Reshape into sum(nTri) x nPts
            vertPos     = reshape(vertPos, nPts, [])';
            % Triangles for element e located in triPos(e):triPos(e+1)-1
            triPos      = [0; cumsum(nTri)] + 1;
            % Make triangulation struct
            tri = struct('Points'          , coords , ...
                         'ConnectivityList', vertPos, ...
                         'triPos'          , triPos , ...
                         'nTri'            , nTri   );
            
        end
        
        %-----------------------------------------------------------------%
        function x = mapCoords(cubature, x, xR)
            % Map cubature points from reference to physical coordinates
            
            G    = cubature.G;
            % Number of triangle vertices
            nPts = 3;
            % Total number of cubature points
            nq   = size(x,1);
            % Total number of triangles
            nTri = sum(cubature.triangulation.nTri);
            % Triangulation points
            ind  = reshape(cubature.triangulation.ConnectivityList', [], 1);
            xt   = cubature.triangulation.Points(ind,:);
            % Compute mappings from reference to physical triangles
            XI = [xR, ones(3,1)];
            X  = zeros(nPts, G.griddim*nTri);
            for dNo = 1:G.griddim
                X(:,dNo:G.griddim:end) = reshape(xt(:,dNo), nPts, []);
            end
            R = XI\X;
            % Map coordinates
            xx = x*R(1:2,:) + R(3,:);
            x  = zeros(nTri*nq, G.griddim);
            for dNo = 1:G.griddim
                xtmp     = xx(:, dNo:G.griddim:end);
                x(:,dNo) = xtmp(:);
            end
            
        end
        
         %-----------------------------------------------------------------%
        function x = mapCoordsBary(cubature, x)
            % Map cubature points from Barycentric to physical coordinates
            
            G    = cubature.G;
            % Number of triangle vertices
            nPts = 3;
            % Total number of cubature points
            nq   = size(x,1);
            % Total number of triangles
            nTri = sum(cubature.triangulation.nTri);
            % Triangulation points
            ind  = reshape(cubature.triangulation.ConnectivityList', [], 1);
            xt   = cubature.triangulation.Points(ind,:);
            
            % Create mapping
            R    = zeros(nPts, G.griddim*nTri);
            for dNo = 1:G.griddim
                R(:,dNo:G.griddim:end) = reshape(xt(:,dNo), nPts, []);
            end
            % Map coordinates
            xx = x*R;
            x  = zeros(nTri*nq, G.griddim);
            for dNo = 1:G.griddim
                xtmp     = xx(:, dNo:G.griddim:end);
                x(:,dNo) = xtmp(:);
            end
            
        end
        
        %-----------------------------------------------------------------%
        function areas = getTriAreas(cubature)
            % Get triangle areas
            
            % Get vector defining two of the edges
            ind = cubature.triangulation.ConnectivityList;
            x   = cubature.triangulation.Points;
            vec1 = x(ind(:,1),:)  - x(ind(:,2),:);
            vec2 = x(ind(:,1),:)  - x(ind(:,3),:);
            % Compute areas
            if cubature.G.griddim == 2
                areas = abs(vec1(:,1).*vec2(:,2) - vec1(:,2).*vec2(:,1))/2;
            else
                areas = sqrt(sum(cross(vec1, vec2, 2).^2, 2))/2;
            end
            
        end
        
        %-----------------------------------------------------------------%
        function cubature = makeCubature(cubature)
            % Make cubature
            
            % Get triangluation of all cells
            cubature.triangulation = cubature.getTriangulation(cubature.G);
            % Number of triangles per cell/face
            nTri = sum(cubature.triangulation.nTri);
            % Cubature points and weights
            [x, w, n, xR] = getTriangleCubaturePointsAndWeights(cubature.prescision);
            % Map to physical coords
            x = cubature.mapCoords(x, xR);
            % Multiply weights by triangle areas
            triNo = reshape(repmat(1:sum(cubature.triangulation.nTri), n, 1), [], 1);
            areas = cubature.getTriAreas();
            w = repmat(w, nTri, 1).*areas(triNo);
            n = cubature.triangulation.nTri.*n;
            cubature = cubature.assignCubaturePointsAndWeights(x,w,n);
            
        end 
            
    end
    
    
    
end

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

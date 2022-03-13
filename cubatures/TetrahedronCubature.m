classdef TetrahedronCubature < Cubature
    
    properties
        triangulation   % Struct containing triangulation information
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function cubature = TetrahedronCubature(G, prescision)
            % Set up tetrahedron cubatureature
            
            % Construction handled by base class
            cubature = cubature@Cubature(G, prescision, 3);
        end
        
        %-----------------------------------------------------------------%
        function tri = getTriangulation(cubature, G) %#ok
            % Triangulate all cells of G
            
            % All faces of all cells
            faces = G.cells.faces(:, 1);
            % All edges of all faces of all cells
            edges = G.faces.edges(mcolon(G.faces.edgePos(faces)     , ...
                                         G.faces.edgePos(faces+1)-1));
            % All nodes of all edges of all faces of all cells
            nodes = G.edges.nodes(mcolon(G.edges.nodePos(edges)     , ...
                                         G.edges.nodePos(edges+1)-1));
            % Number of tetrehedreons belonging to each face of each cell
            nTetf  = diff(G.faces.edgePos);
            nTetf  = nTetf(faces);
            % Number of tetrahedrons per cell
            nTet   = accumarray(rldecode((1:G.cells.num)', diff(G.cells.facePos), 1), nTetf);
            % Coordinates for nodes, faces and cells
            faces  = rldecode(faces, nTetf, 1);
            cells  = rldecode((1:G.cells.num)', nTet, 1);
            xn     = G.nodes.coords(nodes, :);
            xf     = G.faces.centroids(faces,:);
            xc     = G.cells.centroids(cells, :);
            coords = [xn;xf;xc];
            % Number of vertices for a tetrahedron
            nPts     = 4;
            numNodes = numel(nodes);
            numFaces = numel(faces);
            numCells = numel(cells);
            % Index into coordinates for each tetrahedron
            vertPos  = zeros(sum(nTet)*nPts, 1);
            % First two coordinates of each tetrahedron are nodes
            vec         = 1:nPts:sum(nTet)*nPts;
            ii          = mcolon(vec, vec+1);
            vertPos(ii) = 1:numNodes;
            % Third coordinate is face centroid
            ii          = nPts-1:nPts:sum(nTet)*nPts;
            vertPos(ii) = (1:numFaces) + numNodes;
            % Fourth coordinate is cell centroid
            ii          = nPts:nPts:sum(nTet)*nPts;
            vertPos(ii) = (1:numCells) + numNodes + numCells;
            % Reshape into sum(nTet) x nPts
            vertPos = reshape(vertPos, nPts, [])';
            % tetrahedra for element e located in triPos(e):triPos(e+1)-1
            tetPos = [0; cumsum(nTet)] + 1;
            % Make triangulation struct
            tri = struct('Points'          , coords , ...
                         'ConnectivityList', vertPos, ...
                         'triPos'          , tetPos , ...
                         'nTet'            , nTet   );
            
        end
       
        %-----------------------------------------------------------------%
        function x = mapCoords(cubature, x, xR)
            % Map cubature points from Barycentric to physical coordinates
            
            G    = cubature.G;
            % Number of triangle vertices
            nPts = 4;
            % Total number of cubature points
            nq   = size(x,1);
            % Total number of triangles
            nTet = sum(cubature.triangulation.nTet);
            % Triangulation points
            ind  = reshape(cubature.triangulation.ConnectivityList', [], 1);
            xt   = cubature.triangulation.Points(ind,:);
            % Compute mappings from reference to physical triangles
            XI = [xR, ones(4,1)];
            X  = zeros(nPts, G.griddim*nTet);
            for dNo = 1:G.griddim
                X(:,dNo:G.griddim:end) = reshape(xt(:,dNo), nPts, []);
            end
            R = XI\X;
            % Map coordinates
            xx = x*R(1:3,:) + R(4,:);
            x  = zeros(nTet*nq, G.griddim);
            for dNo = 1:G.griddim
                xtmp     = xx(:, dNo:G.griddim:end);
                x(:,dNo) = xtmp(:);
            end
            
        end
        
        %-----------------------------------------------------------------%
        function x = mapCoordsBary(cubature, x)
            % Map cubature points from Barycentric to physical coordinates
            
            G    = cubature.G;
            % Number of tetrahedron vertices
            npts = 4;
            % total number of cubature points
            nq   = size(x,1);
            % Total number of tetrahedrons
            nTet = sum(cubature.triangulation.nTet);
            % Triangulation points
            ind  = reshape(cubature.triangulation.ConnectivityList', [], 1);
            xt   = cubature.triangulation.Points(ind,:);
            % Create mapping
            R    = zeros(npts, G.griddim*nTet);
            for dNo = 1:G.griddim
                R(:,dNo:G.griddim:end) = reshape(xt(:,dNo), npts, []);
            end
            % Map coordinates
            xx = x*R;
            x  = zeros(nTet*nq, G.griddim);
            for dNo = 1:G.griddim
                xtmp     = xx(:, dNo:G.griddim:end);
                x(:,dNo) = xtmp(:);
            end
            
        end
        
        %-----------------------------------------------------------------%
        function volumes = getTetVolumes(cubature)
            % Get tetrahedron volumes
            
            % Get vectors defining three edges
            ind = cubature.triangulation.ConnectivityList;
            x   = cubature.triangulation.Points;
            vec1 = x(ind(:,1),:) - x(ind(:,2),:);
            vec2 = x(ind(:,1),:) - x(ind(:,3),:);
            vec3 = x(ind(:,1),:) - x(ind(:,4),:);
            % Compute volumes
            volumes = abs(dot(cross(vec1, vec2, 2), vec3, 2))/6;
            
        end
        
        %-----------------------------------------------------------------%
        function cubature = makeCubature(cubature)
            % Make cubature
            
            % Get triangulation of all cells
            cubature.triangulation = cubature.getTriangulation(cubature.G);
            % Numbe of tetrahedrons per cell
            nTri = sum(cubature.triangulation.nTet);
            % Cubature points and weights
            [x, w, n, xR] = getTetrahedronCubaturePointsAndWeights(cubature.prescision);
            % Map to physical coords
            x = cubature.mapCoords(x, xR);
            % Multiply weights by tetrahedron volumes
            tetNo     = reshape(repmat(1:sum(cubature.triangulation.nTet), n, 1), [], 1);
            volumes   = cubature.getTetVolumes();
            w = repmat(w, nTri, 1).*volumes(tetNo);
            n = cubature.triangulation.nTet.*n;
            cubature = cubature.assignCubaturePointsAndWeights(x,w,n);
            
        end 
        
        %-----------------------------------------------------------------%
        function plotTriangulation(cubature)
            % Convenience function for plotting triangulation
            
            x  = cubature.triangulation.Points;
            c  = cubature.triangulation.ConnectivityList;
            GT = tetrahedralGrid(x, c);
            
            hold on
            plotGrid(cubature.G)
            plotGrid(GT, 'facec', 'none', 'edgec', 'b');
            hold off
            
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

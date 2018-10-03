classdef TetrahedronCubature < Cubature
    
    properties
        
        parentPos
        tetNo
        triangulation
        
    end
    
    methods
        
        function cub = TetrahedronCubature(G, prescision, internalConn)
            
            cub = cub@Cubature(G, prescision, internalConn);
            cub.triangulation = cub.getTriangulation(G);
            [x, w, n] = cub.makeCubature();
            cub.points = x;
            cub.weights = w;
            cub.numPoints = n;
            
            cub.parentPos = [0; cumsum(cub.triangulation.nTet*cub.numPoints)] + 1;
            
        end
        
        function tri = getTriangulation(cub, G)
            
            faces = G.cells.faces(:, 1);
            edges = G.faces.edges(mcolon(G.faces.edgePos(faces), G.faces.edgePos(faces+1)-1));
            nodes = G.edges.nodes(mcolon(G.edges.nodePos(edges), G.edges.nodePos(edges+1)-1));
            
            nTetf = diff(G.faces.edgePos);
            nTetf = nTetf(faces);
            faces = rldecode(faces, nTetf, 1);
            nTet = accumarray(rldecode((1:G.cells.num)', diff(G.cells.facePos), 1), nTetf);
            parents = rldecode((1:G.cells.num)', nTet, 1);
            
            xn = G.nodes.coords(nodes, :);
            xf = G.faces.centroids(faces,:);
            xc = G.cells.centroids(parents, :);
            coords = [xn;xf;xc];
            
            npts = 4;
            numNodes   = numel(nodes);
            numFaces   = numel(faces);
            numParents = numel(parents);
            
            vertPos = zeros(sum(nTet)*npts, 1);
            
            vec = 1:npts:sum(nTet)*npts;
            ii = mcolon(vec, vec+1);
            vertPos(ii) = 1:numNodes;
            
            ii = npts-1:npts:sum(nTet)*npts;
            vertPos(ii) = (1:numFaces) + numNodes;
            
            ii = npts:npts:sum(nTet)*npts;
            vertPos(ii) = (1:numParents) + numNodes + numParents;
            
            vertPos = reshape(vertPos, npts, [])';
            
            triPos = [0; cumsum(nTet)] + 1;
            
            tri = struct('Points' , coords , ...
                         'ConnectivityList', vertPos, ...
                         'triPos' , triPos , ...
                         'nTet'   , nTet);
            
        end
        
        function plotTriangulation(cub)
            
            x = cub.triangulation.Points;
            c = cub.triangulation.ConnectivityList;
            
            figure;
            hold on
            plotGrid(cub.G)
            GT = tetrahedralGrid(x, c);
            plotGrid(GT, 'facec', 'none', 'edgec', 'b');
        end
        
        function [x, w, n] = getCubaturePointsAndWeights(cub)
            
            % Cubature rules in barycentric coordinates
            
            presc = cub.prescision;
            
            if presc <= 1
                
                w = 1;
                
                x = [0,0,0];
            
            elseif presc == 2
                
                a = 0.138196601125010515179541316563436;
                b = 1 - 3*a;
                
                w = [1;
                     1;
                     1;
                     1]/4;
                 
                x = [b a a a;
                     a b a a;
                     a a b a;
                     a a a b];
                 
            elseif presc == 3
                
                w = [-0.800000000000000
                      0.450000000000000
                      0.450000000000000
                      0.450000000000000
                      0.450000000000000];
                
                x = [0.250000000000000   0.250000000000000   0.250000000000000   0.250000000000000;
                     0.166666666666667   0.500000000000000   0.166666666666667   0.166666666666667;
                     0.500000000000000   0.166666666666667   0.166666666666667   0.166666666666667;
                     0.166666666666667   0.166666666666667   0.166666666666667   0.500000000000000;
                     0.166666666666667   0.166666666666667   0.500000000000000   0.166666666666667];
                 
            elseif presc == 4
                
                
                w = [-0.078933333333333;
                      0.045733333333333;
                      0.045733333333333;
                      0.045733333333333;
                      0.045733333333333;
                      0.149333333333333;
                      0.149333333333333;
                      0.149333333333333;
                      0.149333333333333;
                      0.149333333333333;
                      0.149333333333333];
                
                x = [0.250000000000000   0.250000000000000   0.250000000000000   0.250000000000000;
                     0.071428571428572   0.785714285714286   0.071428571428571   0.071428571428571;
                     0.785714285714286   0.071428571428571   0.071428571428571   0.071428571428571;
                     0.071428571428572   0.071428571428571   0.071428571428571   0.785714285714286;
                     0.071428571428572   0.071428571428571   0.785714285714286   0.071428571428571;
                     0.100596423833201   0.100596423833201   0.399403576166799   0.399403576166799;
                     0.100596423833201   0.399403576166799   0.100596423833201   0.399403576166799;
                     0.100596423833201   0.399403576166799   0.399403576166799   0.100596423833201;
                     0.399403576166799   0.399403576166799   0.100596423833201   0.100596423833201;
                     0.399403576166799   0.100596423833201   0.399403576166799   0.100596423833201;
                     0.399403576166799   0.100596423833201   0.100596423833201   0.399403576166799];

            
            else
                
                error('Prescision not supported')
                
            end
            
            n = numel(w);
            
        end
        
        function x = mapCoords(cub, x)
           
            G = cub.G;
            npts = 4;
            nq = size(x,1);
            nTet = sum(cub.triangulation.nTet);
            ind = reshape(cub.triangulation.ConnectivityList', [], 1);
            xt   = cub.triangulation.Points(ind,:);
            R = zeros(npts, G.griddim*nTet);
            for dNo = 1:G.griddim
                R(:,dNo:G.griddim:end) = reshape(xt(:,dNo), npts, []);
            end
            
            xx = x*R;
            x = zeros(nTet*nq, G.griddim);
            for dNo = 1:G.griddim
                xtmp = xx(:, dNo:G.griddim:end);
                x(:,dNo) = xtmp(:);
            end
            
        end
        
        function volumes = getTetVolumes(cub)
            
            ind = cub.triangulation.ConnectivityList;
            x   = cub.triangulation.Points;
            vec1 = x(ind(:,1),:)  - x(ind(:,2),:);
            vec2 = x(ind(:,1),:)  - x(ind(:,3),:);
            vec3 = x(ind(:,1),:)  - x(ind(:,4),:);
            
            volumes = abs(dot(cross(vec1, vec2, 2), vec3, 2))/6;
            
        end
        
        function [x, w, n, tetNo] = makeCubature(cub)
            
            nTri = sum(cub.triangulation.nTet);
            [x, w, n] = cub.getCubaturePointsAndWeights();
            x = cub.mapCoords(x);
            
            tetNo = reshape(repmat(1:sum(cub.triangulation.nTet), n, 1), [], 1);
            volumes = cub.getTetVolumes();
            
            w = repmat(w, nTri, 1).*volumes(tetNo);
            
        end 
            
    end
    
    
    
end
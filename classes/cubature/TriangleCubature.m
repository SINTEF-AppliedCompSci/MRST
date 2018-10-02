classdef TriangleCubature < Cubature
    
    properties
        
        parentPos
        triNo
        triangulation
        
    end
    
    methods
        
        function cub = TriangleCubature(G, prescision, internalConn)
            
            cub = cub@Cubature(G, prescision, internalConn);
            cub.triangulation = cub.getTriangulation(G);
            [x, w, n] = cub.makeCubature();
            cub.points = x;
            cub.weights = w;
            cub.numPoints = n;
            cub.dim = 2;
            
            cub.parentPos = [0; cumsum(cub.triangulation.nTri*cub.numPoints)] + 1;
            
        end
        
        function tri = getTriangulation(cub, G)
            
            if G.griddim == 2
                faces = G.cells.faces(:, 1);
                nodes = G.faces.nodes(mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces+1)-1));
                nTri  = diff(G.cells.facePos);
                parents = rldecode((1:G.cells.num)', nTri, 1);
                xc    = G.cells.centroids(parents, :);
            else
                edges = G.faces.edges;
                nodes = G.edges.nodes(mcolon(G.edges.nodePos(edges), G.edges.nodePos(edges+1)-1));
                nTri  = diff(G.faces.edgePos);
                parents = rldecode((1:G.faces.num)', nTri, 1);
                xc    = G.faces.centroids(parents, :);
            end

            xn    = G.nodes.coords(nodes, :);
            coords = [xn; xc];
            
            npts = 3;
            numNodes = numel(nodes);
            numParents = numel(parents);
            
            vertPos = zeros(sum(nTri)*npts, 1);
            
            vec = 1:npts:sum(nTri)*npts;
            ii = mcolon(vec, vec+1);
            vertPos(ii) = 1:numNodes;
            
            ii = npts:npts:sum(nTri)*npts;
            vertPos(ii) = (1:numParents) + numNodes;
            
            vertPos = reshape(vertPos, npts, [])';
            
            triPos = [0; cumsum(nTri)] + 1;
            
            tri = struct('Points' , coords , ...
                         'ConnectivityList', vertPos, ...
                         'triPos' , triPos , ...
                         'nTri'   , nTri);
            
        end
        
        function plotTriangulation(cub)
            
            x = cub.triangulation.Points;
            c = cub.triangulation.ConnectivityList;
            
            figure;
            hold on
            plotGrid(cub.G)
            triplot(c, x(:,1), x(:,2));
            
        end
        
        function [x, w, n] = getCubaturePointsAndWeights(cub)
            
            % Cubature rules in barycentric coordinates
            
            presc = cub.prescision;
            
            if presc <= 1
                
                w = 1;
                x = [0,0,0];
            
            elseif presc == 2
                
                w = [1;
                     1;
                     1]/3;
                
                x = [0,1,1;
                     1,0,1;
                     1,1,0]/2;
                
            elseif presc == 3
                
                w = [-0.56250000000000000000;
                      0.52083333333333333333;
                      0.52083333333333333333;
                      0.52083333333333333333];
                  
                  
                x = [0.333333333333333   0.333333333333333   0.333333333333333;
                     0.200000000000000   0.600000000000000   0.200000000000000;
                     0.200000000000000   0.200000000000000   0.600000000000000;
                     0.600000000000000   0.200000000000000   0.200000000000000];
                
            elseif presc == 4

                w = [0.22338158967801;
                     0.22338158967801;
                     0.22338158967801;
                     0.10995174365532;
                     0.10995174365532;
                     0.10995174365532];
                 
                x = [0.108103018168060   0.445948490915970   0.445948490915970;
                     0.445948490915960   0.445948490915970   0.108103018168070;
                     0.445948490915960   0.108103018168070   0.445948490915970;
                     0.816847572980460   0.091576213509770   0.091576213509770;
                     0.091576213509770   0.091576213509770   0.816847572980460;
                     0.091576213509770   0.816847572980460   0.091576213509770];
                 
            elseif presc == 5
                
                w = [0.225000000000000;
                     0.125939180544827;
                     0.125939180544827;
                     0.125939180544827;
                     0.132394152788506;
                     0.132394152788506;
                     0.132394152788506];
                
                x = [0.333333333333333   0.333333333333333   0.333333333333333;
                     0.101286507323456   0.797426985353087   0.101286507323456;
                     0.101286507323456   0.101286507323456   0.797426985353087;
                     0.797426985353087   0.101286507323456   0.101286507323456;
                     0.470142064105115   0.059715871789770   0.470142064105115;
                     0.470142064105115   0.470142064105115   0.059715871789770;
                     0.059715871789770   0.470142064105115   0.470142064105115];
                 
            elseif presc == 6
                
                w =  [0.205950504760887;
                      0.205950504760887;
                      0.205950504760887;
                      0.063691414286223;
                      0.063691414286223;
                      0.063691414286223;
                      0.063691414286223;
                      0.063691414286223;
                      0.063691414286223];
                  
                x = [0.437525248383384   0.124949503233232   0.437525248383384;
                     0.437525248383384   0.437525248383384   0.124949503233232;
                     0.124949503233232   0.437525248383384   0.437525248383384;
                     0.037477420750088   0.797112651860071   0.165409927389841;
                     0.165409927389841   0.797112651860071   0.037477420750088;
                     0.037477420750088   0.165409927389841   0.797112651860071;
                     0.797112651860071   0.165409927389841   0.037477420750088;
                     0.165409927389841   0.037477420750088   0.797112651860071;
                     0.797112651860071   0.037477420750088   0.165409927389841];
                 
            else
                
                error('Prescision not supported')
                
            end
            
            n = numel(w);
            
        end
        
        function x = mapCoords(cub, x)
           
            G = cub.G;
            npts = 3;
            nq = size(x,1);
            nTri = sum(cub.triangulation.nTri);
            ind = reshape(cub.triangulation.ConnectivityList', [], 1);
            xt   = cub.triangulation.Points(ind,:);
            R = zeros(npts, G.griddim*nTri);
            for dNo = 1:G.griddim
                R(:,dNo:G.griddim:end) = reshape(xt(:,dNo), npts, []);
            end
            
            xx = x*R;
            x = zeros(nTri*nq, G.griddim);
            for dNo = 1:G.griddim
                xtmp = xx(:, dNo:G.griddim:end);
                x(:,dNo) = xtmp(:);
            end
            
        end
        
        function areas = getTriAreas(cub)
            
            ind = cub.triangulation.ConnectivityList;
            x   = cub.triangulation.Points;
            vec1 = x(ind(:,1),:)  - x(ind(:,2),:);
            vec2 = x(ind(:,1),:)  - x(ind(:,3),:);
            
            if cub.G.griddim == 2
                areas = abs(vec1(:,1).*vec2(:,2) - vec1(:,2).*vec2(:,1))/2;
            else
                areas = sqrt(sum(cross(vec1, vec2, 2).^2, 2))/2;
            end
            
        end
        
        function [x, w, n, triNo] = makeCubature(cub)
            
            nTri = sum(cub.triangulation.nTri);
            [x, w, n] = cub.getCubaturePointsAndWeights();
            x = cub.mapCoords(x);
            
            triNo = reshape(repmat(1:sum(cub.triangulation.nTri), n, 1), [], 1);
            areas = cub.getTriAreas();
            
            w = repmat(w, nTri, 1).*areas(triNo);
            
        end 
            
    end
    
    
    
end
classdef UnstructCubature < Cubature
    
    properties
        
        parentPos
        
    end
    
    methods
        
        function cub = UnstructCubature(G, prescision)
            
            cub = cub@Cubature(G, prescision);
            [x, w, n] = cub.makeCubature();
            cub.points = x;
            cub.weights = w;
            cub.numPoints = n;
            
            
            cub.parentPos = (0:cub.numPoints:G.cells.num*cub.numPoints)' + 1;
%             cub.parentPos = [0; cumsum(cub.triangulation.nTri*cub.numPoints)] + 1;
            
        end
        
        function [x, w, n] = getCubaturePointsAndWeights(cub)
            
            presc = cub.prescision;
            
            if presc <= 1
                
                w = 1;
                x = zeros(1, cub.G.griddim);
             
            else
                
                error('Prescision not supported')
                
            end
            
            n = numel(w);
            
        end
        
        function x = mapCoords(cub, x)
           
            G = cub.G;
            nq = size(x,1);
            
            xc = rldecode(G.cells.centroids, nq, 1);
            x = x + xc;
            
        end
        
        function [x, w, n, cNo] = makeCubature(cub)
            
            G = cub.G;
            [x, w, n] = cub.getCubaturePointsAndWeights();
            x = cub.mapCoords(x);
            
            cNo = rldecode((1:G.cells.num)', n, 1);
            volumes = G.cells.volumes(cNo);
            
            w = repmat(w, G.cells.num*n, 1).*volumes;
            
        end 
        
    end
    
end
classdef CellTriangleCubature < TriangleCubature
    
    methods
        
        function cub = CellTriangleCubature(G, prescision)
            
            cub = cub@TriangleCubature(G, prescision);
            cub.triangulation = cub.getTriangulation(G);
            [x, w, n] = cub.makeCubature();
            cub.points = x;
            cub.weights = w;
            cub.numPoints = n;
            
            cub.parentPos = [0; cumsum(cub.triangulation.nTri*cub.numPoints)] + 1;
            
        end
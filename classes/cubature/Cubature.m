classdef Cubature
    
    properties
        
        G
        prescision
        points
        weights
        numPoints
        
    end
    
    methods
        
        function cub = Cubature(G, prescision)
            
            cub.G = G;
            cub.prescision = prescision;
            cub.points  = [];
            cub.weights = [];
            
        end
            
    end
    
end
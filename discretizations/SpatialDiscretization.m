classdef SpatialDiscretization
        
    properties
        G
        operators
        internalConn
        N
        suffix = '';
    end
    
    methods 
        function disc = SpatialDiscretization(model)
            G                 = model.G;
            disc.G            = G;
            N                 = G.faces.neighbors;
            disc.internalConn = all(N ~= 0, 2);
            disc.N            = N(disc.internalConn,:);
            disc.operators    = model.operators;
        end
        
    end
    
end
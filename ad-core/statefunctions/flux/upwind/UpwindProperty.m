classdef UpwindProperty
    % Base class which adds the faceUpstream function to derived classes
    properties
        UpwindDiscretization % Class instance
    end
    
    methods
        function up = UpwindProperty(upstream)
            up.UpwindDiscretization = upstream;
        end
        
        function v = faceUpstream(prop, model, state, flag, cellvalue)
            v = prop.UpwindDiscretization.faceUpstream(model, state, flag, cellvalue);
        end
    end
end
classdef UpwindProperty
    properties
        UpstreamDiscretization
    end
    
    methods
        function up = UpwindProperty(upstream)
            up.UpstreamDiscretization = upstream;
        end
        
        function v = faceUpstream(prop, state, flag, cellvalue)
            v = prop.UpstreamDiscretization.faceUpstream(state, flag, cellvalue);
        end
    end
end
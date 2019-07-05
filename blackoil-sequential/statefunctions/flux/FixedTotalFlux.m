classdef FixedTotalFlux < StateFunction
    properties
    end
    
    methods
        function gp = FixedTotalFlux(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn({'flux'}, 'state');
        end
        function vT = evaluateOnDomain(prop, model, state)
            vT = sum(state.flux(model.operators.internalConn, :), 2);
        end
    end
end
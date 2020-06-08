classdef BiotGradP < StateFunction
    
    methods
        function gp = BiotGradP(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn('biotgradp', 'state');
            gp.label = '\alpha \nabla p';
        end
        
        function biotgradp = evaluateOnDomain(prop, model, state)
            biotgradp = model.getProp(state, 'biotgradp');
        end
    end
    
end

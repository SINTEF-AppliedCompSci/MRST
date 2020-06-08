classdef Stress < StateFunction
    
    methods
        function gp = Stress(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn('displacement', 'state');
            gp.label = 'stress';
        end
        
        function stress = evaluateOnDomain(prop, model, state)
            op = model.operators.global_stress;
            u = model.getProp(state, 'displacement');
            error('not yet set up for mpsa');
            stress = op*u;
        end
    end
    
end

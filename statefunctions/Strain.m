classdef Strain < StateFunction
    
    methods
        function gp = Strain(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn('displacement', 'state');
            gp.label = 'strain';
        end
        
        function strain = evaluateOnDomain(prop, model, state)
            op = model.operators.global_strain;
            u = model.getProp(state, 'displacement');
            error('not yet implemented for mpsa')
            strain = op*u;
        end
    end
    
end

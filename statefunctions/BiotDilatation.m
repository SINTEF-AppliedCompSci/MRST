classdef BiotDilatation < StateFunction
    
    methods
        function gp = BiotDilatation(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn('dilatation', 'state');
            gp.label = '\nabla\cdot u';
        end
        
        function divu = evaluateOnDomain(prop, model, state)
            divu = model.getProp(state, 'dilatation');
        end
        
    end
    
end

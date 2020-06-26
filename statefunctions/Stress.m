classdef Stress < StateFunction
    
    methods
        function gp = Stress(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn('FaceNodeDisplacement');
            gp = gp.dependsOn('displacement', 'state');
            gp.label = 'stress';
        end
        
        function stress = evaluateOnDomain(prop, model, state)
            stressop = model.operators.stressop;
            [unf, uc] = model.getProps(state, 'FaceNodeDisplacement', 'displacement');
            stress = stressop(unf, uc);
        end
    end
    
end

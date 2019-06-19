classdef SurfactantPhasePressures < StateFunction
    properties
    end
    
    methods
        function gp = SurfactantPhasePressures(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn({'CapillaryPressure'});
            gp = gp.dependsOn({'pressure'}, 'state');
            gp = gp.dependsOn({'surfactant'}, 'state');
        end
        
        function p_phase = evaluateOnDomain(prop, model, state)
            fluid = model.fluid;
            cs = model.getProps(state, 'surfactant');
            p = model.getProps(state, 'Pressure');
            pc = prop.getEvaluatedDependencies(state, 'CapillaryPressure');
            pc{1} = pc{1}.*fluid.ift(cs)/fluid.ift(0);
            pc_phase = pc;
        end
    end
end
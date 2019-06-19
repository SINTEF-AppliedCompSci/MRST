classdef SurfactantCapillaryPressure < StateFunction
    properties
    end
    
    methods
        function gp = SurfactantCapillaryPressure(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn({'CapillaryPressure'});
            gp = gp.dependsOn({'pressure'}, 'state');
            gp = gp.dependsOn({'surfactant'}, 'state');
        end
        
        function p_phase = evaluateOnDomain(prop, model, state)
            fluid = model.fluid;
            c = model.getProps(state, 'surfactant');
            p = model.getProps(state, 'Pressure');
            pcow = prop.evaluateFunctionOnDomainWithArguments(fluid.pcOW, sW);
            pcow = pcow.*fluid.ift(c)/fluid.ift(0);
            pc_phase = pc;
        end
    end
end
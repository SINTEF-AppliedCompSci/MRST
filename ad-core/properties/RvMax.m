classdef RvMax < GridProperty
    properties
    end
    
    methods
        function rsSat = evaluateOnGrid(prop, model, state)
            p = model.getProp(state, 'pressure');
            pc = state.FlowProps.CapillaryPressure{model.water + model.oil + model.gas};
            if model.disgas
                f = model.fluid;
                rsSat = prop.evaluateFunctionOnGrid(f.rsSat, p + pc);
            else
                rsSat = 0*p;
            end
        end
    end
end
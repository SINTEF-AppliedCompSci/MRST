classdef RvMax < GridProperty
    properties
    end
    
    methods
        function rsSat = evaluateOnGrid(prop, model, state)
            p = model.getProp(state, 'pressure');
            pc = state.FlowProps.CapillaryPressure{model.water + model.oil + model.gas};
            if model.disgas
                f = model.fluid;
                pg = p;
                if ~isempty(pc)
                    pg = pg + pc;
                end
                rsSat = prop.evaluateFunctionOnGrid(f.rvSat, pg);
            else
                rsSat = 0*p;
            end
        end
    end
end
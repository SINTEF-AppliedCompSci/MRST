classdef MultipliedPoreVolume < GridProperty
    properties
    end
    
    methods
        function pv = evaluateOnGrid(prop, model, state)
            f = model.fluid;
            pv = model.operators.pv;
            if isfield(f, 'pvMultR')
                p = model.getProp(state, 'pressure');
                pvMult = prop.evaluateFunctionOnGrid(f.pvMultR, p);
                pv = pv.*pvMult;
            end
        end
    end
end
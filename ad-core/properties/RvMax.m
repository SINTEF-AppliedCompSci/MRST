classdef RvMax < GridProperty
    properties
        rvReduction = 0;
    end
    
    methods
        function rvSat = evaluateOnDomain(prop, model, state)
            p = model.getProp(state, 'pressure');
            pc = state.FlowProps.CapillaryPressure{model.water + model.oil + model.gas};
            if model.disgas
                f = model.fluid;
                pg = p;
                if ~isempty(pc)
                    pg = pg + pc;
                end
                rvSat = prop.evaluateFunctionOnGrid(f.rvSat, pg);
                if prop.rvReduction > 0 && isfield(state, 'sMax')
                    [sOMax, sO] = model.getProps(state, 'somax', 'so');
                    factor = (sOMax./sO);
                    factor(value(sO) == 0) = 1;
                    rvSat = rvSat.*factor.^prop.rvReduction;
                end
            else
                rvSat = 0*p;
            end
        end
    end
end
classdef RvMax < GridProperty
    properties
        rvReduction = 0;
    end
    
    methods
        function gp = RvMax(varargin)
            gp@GridProperty(varargin{:});
        end
        
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
                    sOMax = max(sOMax, sO);
                    factor = (sO + 1e-4)./(sOMax + 1e-4);
                    rvSat = rvSat.*(factor.^prop.rvReduction);
                end
            else
                rvSat = 0*p;
            end
        end
    end
end
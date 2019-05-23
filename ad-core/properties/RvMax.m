classdef RvMax < StateFunction
    properties
        rvReduction = 0;
    end
    
    methods
        function gp = RvMax(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn('PhasePressures');
        end
        
        function rvSat = evaluateOnDomain(prop, model, state)
            p_phase = prop.getEvaluatedDependencies(state, 'PhasePressures');
            pg = p_phase{model.water + model.oil + model.gas};
            if model.vapoil
                f = model.fluid;
                rvSat = prop.evaluateFunctionOnGrid(f.rvSat, pg);
                if prop.rvReduction > 0 && isfield(state, 'sMax')
                    [sOMax, sO] = model.getProps(state, 'somax', 'so');
                    sOMax = max(sOMax, sO);
                    factor = (sO + 1e-4)./(sOMax + 1e-4);
                    rvSat = rvSat.*(factor.^prop.rvReduction);
                end
            else
                rvSat = 0*pg;
            end
        end
    end
end
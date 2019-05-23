classdef CompositionalDensity < StateFunction
    properties
    end
    
    methods
        function gp = CompositionalDensity(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PhasePressures', 'PhaseCompressibilityFactors', 'ComponentPhaseMoleFractions'});
        end

        function rho = evaluateOnDomain(prop, model, state)
            [p, T] = model.getProps(state, 'pressure', 'temperature');
            [p_phase, Z, mf] = prop.getEvaluatedDependencies(state, ...
                'PhasePressures', 'PhaseCompressibilityFactors', 'ComponentPhaseMoleFractions');
            hasWater = model.water;
            
            L_ix = 1+model.water;
            V_ix = L_ix + 1;
            
            x = mf((1+model.water):end, L_ix);
            y = mf((1+model.water):end, V_ix);

            rhoL = model.EOSModel.PropertyModel.computeDensity(p, x, Z{L_ix}, T, true);
            rhoV = model.EOSModel.PropertyModel.computeDensity(p, y, Z{V_ix}, T, false);
            
            if hasWater
                f = model.fluid;
                bW = prop.evaluateFunctionOnDomainWithArguments(f.bW, p_phase{1});
                rhoW = f.rhoWS.*bW;
                rho = {rhoW, rhoL, rhoV};
            else
                rho = {rhoL, rhoV};
            end
        end
    end
end
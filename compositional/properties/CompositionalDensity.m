classdef CompositionalDensity < StateFunction
    % Density where the liquid and vapor phases are predicted by a EOS
    properties
        useCompactEvaluation = true;
    end
    
    methods
        function gp = CompositionalDensity(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PhasePressures', 'PhaseCompressibilityFactors', 'ComponentPhaseMoleFractions'});
            gp = gp.dependsOn({'pressure', 'T'}, 'state');
            gp.label = '\rho_\alpha';
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
            pm = model.EOSModel.PropertyModel;
            rhoL = pm.computeDensity(p, x, Z{L_ix}, T, true);
            if prop.useCompactEvaluation
                [~, ~, twoPhase] = model.getFlag(state);
                if all(twoPhase)
                    rhoV = pm.computeDensity(p, y, Z{V_ix}, T, false);
                else
                    rhoV = rhoL;
                    if any(twoPhase)
                        y2ph = cellfun(@(x) x(twoPhase), y, 'UniformOutput', false);
                        rhoV(twoPhase) = pm.computeDensity(p(twoPhase), y2ph, Z{V_ix}(twoPhase), T(twoPhase), false);
                    end
                end
            else
                rhoV = pm.computeDensity(p, y, Z{V_ix}, T, false);
            end
            
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

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

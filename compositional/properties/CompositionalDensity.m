classdef CompositionalDensity < StateFunction
    % Density where the liquid and vapor phases are predicted by a EOS
    properties
        useCompactEvaluation = true;
    end
    
    methods
        function gp = CompositionalDensity(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PhasePressures',...
                               'PhaseCompressibilityFactors',...
                               'ComponentPhaseMoleFractions'});
            gp = gp.dependsOn({'pressure', 'T'}, 'state');
            gp.label = '\rho_\alpha';
        end

        function rho = evaluateOnDomain(prop, model, state)
            [p, T] = model.getProps(state, 'pressure', 'temperature');
            [p_phase, Z, mf] = prop.getEvaluatedDependencies(state, ...
                'PhasePressures', 'PhaseCompressibilityFactors', 'ComponentPhaseMoleFractions');
            phases = model.getPhaseNames();
            nph = numel(phases);
            L_ix = model.getLiquidIndex();
            V_ix = model.getVaporIndex();
            isEoS = model.getEoSComponentMask();
            
            x = mf(isEoS, L_ix);
            y = mf(isEoS, V_ix);
            eos = model.EOSModel;
            pm = eos.PropertyModel;
            rhoL = pm.computeDensity(eos, p, x, Z{L_ix}, T, true);
            if prop.useCompactEvaluation
                [~, ~, twoPhase] = model.getFlag(state);
                if all(twoPhase)
                    rhoV = pm.computeDensity(eos, p, y, Z{V_ix}, T, false);
                else
                    rhoV = rhoL;
                    if any(twoPhase)
                        y2ph = cellfun(@(x) x(twoPhase), y, 'UniformOutput', false);
                        rhoV(twoPhase) = pm.computeDensity(eos, p(twoPhase), y2ph, Z{V_ix}(twoPhase), T(twoPhase), false);
                    end
                end
            else
                rhoV = pm.computeDensity(eos, p, y, Z{V_ix}, T, false);
            end
            rho = cell(1, nph);
            rho{L_ix} = rhoL;
            rho{V_ix} = rhoV;
            for i = 1:nph
                if i == L_ix || i == V_ix
                    continue
                end
                sn = phases(i);
                b = prop.evaluateFluid(model, ['b', sn], p_phase{i});
                rhoS = model.fluid.(['rho', sn, 'S']);
                rho{i} = rhoS.*b;
            end
        end
    end
end

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

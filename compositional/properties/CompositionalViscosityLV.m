classdef CompositionalViscosityLV < Viscosity
    properties
        useCompactEvaluation = true;
    end
    
    methods
        function gp = CompositionalViscosityLV(model, varargin)
            gp@Viscosity(model, varargin{:});
            gp = gp.dependsOn({'PhasePressures', 'PhaseCompressibilityFactors', 'ComponentPhaseMoleFractions'});
            gp = gp.dependsOn({'pressure', 'T'}, 'state');
            gp.label = '\mu_\alpha';
        end
        
        function mu = evaluateOnDomain(prop, model, state)
            ph_names = model.getPhaseNames();
            nph = numel(ph_names);
            mu = cell(1, nph);
            
            [p, T] = model.getProps(state, 'pressure', 'T');
            [Z, mf] = prop.getEvaluatedDependencies(state, 'PhaseCompressibilityFactors',...
                                                           'ComponentPhaseMoleFractions');
            L_ix = model.getLiquidIndex();
            V_ix = model.getVaporIndex();
            isEoS = model.getEoSComponentMask();

            x = mf(isEoS, L_ix);
            y = mf(isEoS, V_ix);

            eos = model.EOSModel;
            pm = eos.PropertyModel;
            mu{L_ix} = pm.computeViscosity(eos, p, x, Z{L_ix}, T, true);
            
            twoPhase = model.getTwoPhaseFlag(state);
            if prop.useCompactEvaluation && ~all(twoPhase)
                muV = mu{L_ix};
                if any(twoPhase)
                    if iscell(y)
                        y = cellfun(@(x) x(twoPhase), y, 'uniformoutput', false);
                    else
                        y = y(twoPhase, :);
                    end
                    muV(twoPhase) = pm.computeViscosity(eos, p(twoPhase), y, Z{V_ix}(twoPhase), T(twoPhase), false);
                end
            else
                muV = pm.computeViscosity(eos, p, y, Z{V_ix}, T, false);
            end
            mu{V_ix} = muV;
            % Deal with non-flash components
            for i = 1:nph
                if i == L_ix || i == V_ix
                    continue
                end
                mu{i} = prop.evaluatePhaseViscosity(model, state, ph_names(i), p);
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

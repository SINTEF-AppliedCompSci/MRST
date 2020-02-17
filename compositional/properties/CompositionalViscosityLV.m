classdef CompositionalViscosityLV < StateFunction
    properties
        useCompactEvaluation = true;
    end
    
    methods
        function gp = CompositionalViscosityLV(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PhasePressures', 'PhaseCompressibilityFactors', 'ComponentPhaseMoleFractions'});
            gp = gp.dependsOn({'pressure', 'temperature'}, 'state');
            gp.label = '\mu_\alpha';
        end
        
        function mu = evaluateOnDomain(prop, model, state)
            [act, phInd] = model.getActivePhases();
            nph = sum(act);
            mu = cell(1, nph);
            
            [p, T] = model.getProps(state, 'pressure', 'T');
            [Z, mf] = prop.getEvaluatedDependencies(state, 'PhaseCompressibilityFactors',...
                                                           'ComponentPhaseMoleFractions');
            oix = phInd == 2;
            gix = phInd == 3;
            wat = model.water;
            x = mf(1:end-wat, oix);
            y = mf(1:end-wat, gix);

            if model.water
                p_phase = prop.getEvaluatedDependencies(state, 'PhasePressures');
                f = model.fluid;
                wix = phInd == 1;
                pw = p_phase{wix};
                mu{wix} = prop.evaluateFunctionOnDomainWithArguments(f.muW, pw);
            end
            eos = model.EOSModel;
            pm = eos.PropertyModel;
            mu{oix} = pm.computeViscosity(eos, p, x, Z{oix}, T, true);
            
            twoPhase = model.getTwoPhaseFlag(state);
            if prop.useCompactEvaluation && ~all(twoPhase)
                muV = mu{oix};
                if any(twoPhase)
                    if iscell(y)
                        y = cellfun(@(x) x(twoPhase), y, 'uniformoutput', false);
                    else
                        y = y(twoPhase, :);
                    end
                    muV(twoPhase) = pm.computeViscosity(eos, p(twoPhase), y, Z{gix}(twoPhase), T(twoPhase), false);
                end
            else
                muV = pm.computeViscosity(eos, p, y, Z{gix}, T, false);
            end
            mu{gix} = muV;
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

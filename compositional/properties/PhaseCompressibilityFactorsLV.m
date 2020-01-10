classdef PhaseCompressibilityFactorsLV < StateFunction
    properties
        useCompactEvaluation = true;
    end
    
    methods
        function gp = PhaseCompressibilityFactorsLV(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PhaseMixingCoefficients', 'ComponentPhaseMoleFractions'});
            gp = gp.dependsOn({'pressure'}, 'state');
        end

        function v = evaluateOnDomain(prop, model, state)
            if isfield(state, 'cellJacMap')
                arg = {state.cellJacMap};
            else
                arg = {};
            end
            nph = model.getNumberOfPhases();
            eos = model.EOSModel;
            p = model.getProps(state, 'pressure');
            [mix, mf] = prop.getEvaluatedDependencies(state, 'PhaseMixingCoefficients', 'ComponentPhaseMoleFractions');
            L_ix = 1+model.water;
            V_ix = L_ix + 1;
            
            x = mf((1+model.water):end, L_ix);
            y = mf((1+model.water):end, V_ix);
            
            L_mix = mix{L_ix};
            V_mix = mix{V_ix};
            s = getSampleAD(p, x{:}, y{:});
            
            Z_L = eos.computeCompressibilityZ(p, x, L_mix.A, L_mix.B, L_mix.Si, L_mix.Bi, true);
            Z_L = model.AutoDiffBackend.convertToAD(Z_L, s);
            Z_L = eos.setZDerivatives(Z_L, L_mix.A, L_mix.B, arg{:});
            twoPhase = model.getTwoPhaseFlag(state);
            if prop.useCompactEvaluation && ~all(twoPhase) 
                Z_V = Z_L;
                if any(twoPhase)
                    p = p(twoPhase);
                    if iscell(y)
                        y = cellfun(@(x) x(twoPhase), y, 'UniformOutput', false);
                    else
                        y = y(twoPhase, :);
                    end
                    if iscell(V_mix.Si)
                        Si = cellfun(@(x) x(twoPhase), V_mix.Si, 'UniformOutput', false);
                        Bi = cellfun(@(x) x(twoPhase), V_mix.Bi, 'UniformOutput', false);
                    else
                        Si = V_mix.Si(twoPhase, :);
                        Bi = V_mix.Bi(twoPhase, :);
                    end
                    Z_V(twoPhase) = eos.computeCompressibilityZ(p, y, V_mix.A(twoPhase), V_mix.B(twoPhase), Si, Bi, false);
                    Z_V(twoPhase) = eos.setZDerivatives(Z_V(twoPhase), V_mix.A(twoPhase), V_mix.B(twoPhase));
                end
            else
                Z_V = eos.computeCompressibilityZ(p, y, V_mix.A, V_mix.B, V_mix.Si, V_mix.Bi, false);
                Z_V = model.AutoDiffBackend.convertToAD(Z_V, s);
                Z_V = eos.setZDerivatives(Z_V, V_mix.A, V_mix.B, arg{:});
            end
            
            v = cell(1, nph);
            v{L_ix} = Z_L;
            v{V_ix} = Z_V;
            if model.water
                v{1} = ones(numelValue(Z_L), 1);
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

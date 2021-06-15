classdef PhaseMixingCoefficientsLV < StateFunction
    % Mixing coefficients for liquid-vapor system
    properties
        useCompactEvaluation = true;
    end
    
    methods
        function gp = PhaseMixingCoefficientsLV(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'pressure', 'T'}, 'state');
            gp = gp.dependsOn('ComponentPhaseMoleFractions');
            gp.label = 'S_i B_i A_\alpha B_\alpha A_{ij}';
        end

        function v = evaluateOnDomain(prop, model, state)
            eos = model.EOSModel;
            [p, T] = model.getProps(state, 'pressure', 'temperature');
            xy = prop.getEvaluatedDependencies(state, 'ComponentPhaseMoleFractions');
            L_ix = model.getLiquidIndex();
            V_ix = model.getVaporIndex();
            isEoS = model.getEoSComponentMask();
            nph = size(xy, 2);
            v = cell(1, nph);
            
            x = xy(isEoS, L_ix);
            y = xy(isEoS, V_ix);
            acf = eos.CompositionalMixture.acentricFactors;
            
            [A_ij, Bi] = eos.getMixingParameters(p, T, acf, iscell(x));
            [Si_L, A_L, B_L] = eos.getPhaseMixCoefficients(x, A_ij, Bi);
            twoPhase = model.getTwoPhaseFlag(state);
            if prop.useCompactEvaluation && iscell(y) && ~all(twoPhase)
                Si_V = Si_L;
                A_V = A_L;
                B_V = B_L;
                if any(twoPhase)
                    A_ij_2ph = cellfun(@(x) x(twoPhase), A_ij, 'UniformOutput', false);
                    Bi_2ph = cellfun(@(x) x(twoPhase), Bi, 'UniformOutput', false);
                    y_2ph = cellfun(@(x) x(twoPhase), y, 'UniformOutput', false);
                    [Si_V_2ph, A_V(twoPhase), B_V(twoPhase)] = eos.getPhaseMixCoefficients(y_2ph, A_ij_2ph, Bi_2ph);
                    for i = 1:numel(Si_V_2ph)
                        Si_V{i}(twoPhase) = Si_V_2ph{i};
                    end
                end
            else
                [Si_V, A_V, B_V] = eos.getPhaseMixCoefficients(y, A_ij, Bi);
            end
            
            v{L_ix} = struct('Si', {Si_L}, 'A', {A_L}, 'B', {B_L}, 'Bi', {Bi}, 'A_ij', {A_ij});
            v{V_ix} = struct('Si', {Si_V}, 'A', {A_V}, 'B', {B_V}, 'Bi', {Bi}, 'A_ij', {A_ij});
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

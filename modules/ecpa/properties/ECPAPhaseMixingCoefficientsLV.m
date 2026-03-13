classdef ECPAPhaseMixingCoefficientsLV < StateFunction
    % Mixing coefficients for liquid-vapor system
    properties
        useCompactEvaluation = true;
    end
    
    methods
        function gp = ECPAPhaseMixingCoefficientsLV(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'T'}, 'state');
            gp = gp.dependsOn('ComponentPhaseMoleFractions');
            gp.label = 'Si A B Ai Tr';
        end

        function v = evaluateOnDomain(prop, model, state)
            eos = model.EOSModel;
            T = model.getProps(state, 'temperature');
            xy = prop.getEvaluatedDependencies(state, 'ComponentPhaseMoleFractions');
            L_ix = model.getLiquidIndex();
            V_ix = model.getVaporIndex();
            isEoS = model.getEoSComponentMask();
            nph = size(xy, 2);
            v = cell(1, nph);
            
            x = xy(isEoS, L_ix);
            y = xy(isEoS, V_ix);
            
            [A_ij, Ai, Tr] = eos.getMixingParameters(T, iscell(x));
            [Si_L, A_L, B_L] = eos.getPhaseMixCoefficients(x, A_ij);
            twoPhase = model.getTwoPhaseFlag(state);
            if prop.useCompactEvaluation && iscell(y) && ~all(twoPhase)
                A_V = A_L;
                B_V = B_L;
                Si_V = Si_L;
                if any(twoPhase)
                    y_2ph = cellfun(@(x) x(twoPhase), y, 'UniformOutput', false);
                    A_ij = cellfun(@(x) x(twoPhase), A_ij, 'UniformOutput', false);
                    [Si, A_V(twoPhase), B_V(twoPhase)] = eos.getPhaseMixCoefficients(y_2ph, A_ij);
                    for i = 1:numel(Si_V)
                        Si_V{i}(twoPhase) = Si{i};
                    end
                end
            else
                [Si_V, A_V, B_V] = eos.getPhaseMixCoefficients(y, A_ij);
            end
            
            v{L_ix} = struct('Si', {Si_L}, 'A', {A_L}, 'B', {B_L}, 'Ai', {Ai}, 'Tr', {Tr});
            v{V_ix} = struct('Si', {Si_V}, 'A', {A_V}, 'B', {B_V}, 'Ai', {Ai}, 'Tr', {Tr});
        end
    end
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

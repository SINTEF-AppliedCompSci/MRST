classdef ECPAPhaseCompressibilityFactorsLV < StateFunction
    % Compressibility factors for liquid-vapor system (multiplier to ideal
    % gas law)
    properties
        useCompactEvaluation = true;
    end
    
    methods
        function gp = ECPAPhaseCompressibilityFactorsLV(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PhaseMixingCoefficients', 'ComponentPhaseMoleFractions'});
            gp = gp.dependsOn({'pressure', 'T'}, 'state');
            gp.label = 'Z_\alpha';
        end
        
        function v = evaluateOnDomain(prop, model, state)
            if isfield(state, 'cellJacMap')
                arg = {state.cellJacMap};
            else
                arg = {};
            end
            nph = model.getNumberOfPhases();
            eos = model.EOSModel;
            [p, T] = model.getProps(state, 'pressure', 'T');
            [mix, mf] = prop.getEvaluatedDependencies(state, 'PhaseMixingCoefficients', 'ComponentPhaseMoleFractions');
            L_ix = model.getLiquidIndex();
            V_ix = model.getVaporIndex();
            isEoS = model.getEoSComponentMask();
            
            x = mf(isEoS, L_ix);
            y = mf(isEoS, V_ix);
            
            L_mix = mix{L_ix};
            V_mix = mix{V_ix};
            
            [~, v_L, XA_L, XC_L] = eos.computeCompressibilityZ(p, T, x, L_mix.A, L_mix.B, L_mix.Ai, L_mix.Si, L_mix.Tr, true);
            [v_L, XA_L, XC_L] = eos.setvDerivatives(T, p, x, L_mix.A, L_mix.B, v_L, XA_L, XC_L, L_mix.Tr, arg{:});
            R = 8.314462618;
            Z_L = p.*v_L./(R .* T);
            twoPhase = model.getTwoPhaseFlag(state);
            if prop.useCompactEvaluation && ~all(twoPhase)
                Z_V = Z_L;
                v_V = v_L;
                XA_V = XA_L;
                XC_V = XC_L;
                if any(twoPhase)
                    isDiagonal = isa(model.AutoDiffBackend, 'DiagonalAutoDiffBackend');
                    if isDiagonal
                        if iscell(y)
                            y = cellfun(@(x) x(twoPhase), y, 'UniformOutput', false);
                        else
                            y = y(twoPhase, :);
                        end
                        if iscell(V_mix.Si)
                            Si = cellfun(@(x) x(twoPhase), V_mix.Si, 'UniformOutput', false);
                            Tr = cellfun(@(x) x(twoPhase), V_mix.Tr, 'UniformOutput', false);
                            Ai = cellfun(@(x) x(twoPhase), V_mix.Ai, 'UniformOutput', false);
                        else
                            Si = V_mix.Si(twoPhase, :);
                            Tr = V_mix.Tr(twoPhase, :);
                            Ai = V_mix.Ai(twoPhase, :);
                        end
                        A = V_mix.A(twoPhase); 
                        B = V_mix.B(twoPhase); 
                        [~, v_V(twoPhase), XA2_V, XC_V(twoPhase)] = eos.computeCompressibilityZ(p(twoPhase), T(twoPhase), y, A, B, Ai, Si, Tr, false);
                        [v_V(twoPhase), XA2_V, XC_V(twoPhase)] = eos.setvDerivatives(T(twoPhase), p(twoPhase), y, A, B, v_V(twoPhase), XA2_V, XC_V(twoPhase), Tr);
  
                        if iscell(XA2_V)
                            for i = 1:numel(XA2_V)
                                XA_V{i}(twoPhase) =  XA2_V{i};
                            end
                        else
                            XA_V(twoPhase,:) = XA2_V;
                        end
                        Z_V = p.*v_V./(R .* T);
                    else
                        p = p(twoPhase);
                        p = prop.dpADI(p, 1, twoPhase);
                        T = T(twoPhase);
                        if iscell(y)
                            y = cellfun(@(x) x(twoPhase), y, 'UniformOutput', false);
                            y = cellfun(@(x) prop.dpADI(x, 1, twoPhase), y, 'UniformOutput', false);
                        else
                            y = y(twoPhase, :);
                        end
                        if iscell(V_mix.acti)
                            Si = cellfun(@(x) x(twoPhase), V_mix.Si, 'UniformOutput', false);
                            Si = cellfun(@(x) prop.dpADI(x, 1, twoPhase), Si, 'UniformOutput', false);
                            Tr = cellfun(@(x) x(twoPhase), V_mix.Tr, 'UniformOutput', false);
                            Ai = cellfun(@(x) x(twoPhase), V_mix.Ai, 'UniformOutput', false);
                        else
                            Si = V_mix.Si(twoPhase, :);
                            Tr = V_mix.Tr(twoPhase, :);
                            Ai = V_mix.Ai(twoPhase, :);
                        end
                        A = V_mix.A(twoPhase); A = prop.dpADI(A, 1, twoPhase);
                        B = V_mix.B(twoPhase); B = prop.dpADI(B, 1, twoPhase);
                        [~, v_V2, XA2_V, XC_V2] = eos.computeCompressibilityZ(p, T, y, A, B, Ai, Si, Tr, false);
                        [v_V2, XA2_V, XC_V2] = eos.setvDerivatives(T, p, y, A, B, v_V2, XA2_V, XC_V2, Tr);
                        Z_V2 = p.*v_V2./(R .* T);
                        num_two = sum(twoPhase);
                        subs = 1:num_two;
                        if iscell(XA2_V)
                            for i = 1:numel(XA2_V)
                                XA_V{i} = prop.replaceADI(XA_V{i}, XA2_V{i}, twoPhase, subs);
                            end
                        else
                            XA_V(twoPhase,:) = XA2_V;
                        end
                        v_V = prop.replaceADI(v_V, v_V2, twoPhase, subs);
                        XC_V = prop.replaceADI(XC_V, XC_V2, twoPhase, subs);
                        Z_V = prop.replaceADI(Z_V, Z_V2, twoPhase, subs);
                    end
                end
            else
                [~, v_V, XA_V, XC_V] = eos.computeCompressibilityZ(p, T, y, V_mix.A, V_mix.B, V_mix.Ai, V_mix.Si, V_mix.Tr, false);
                [v_V, XA_V, XC_V] = eos.setvDerivatives(T, p, y, V_mix.A, V_mix.B, v_V, XA_V, XC_V, V_mix.Tr, arg{:});
                Z_V = p.*v_V./(R * T);
            end
            
            v = cell(1, nph);
            v{L_ix} = struct('Z', {Z_L}, 'v', {v_L},'XA', {XA_L},'XC', {XC_L});
            v{V_ix} = struct('Z', {Z_V}, 'v', {v_V},'XA', {XA_V},'XC', {XC_V});
            if nph > 2
                [v{cellfun(@isempty, v)}] = deal(ones(numelValue(Z_L), 1));
            end
        end
    end
    
    methods (Static, Access=protected)
        function x = dpADI(x, ix, subs)
            if islogical(subs)
                subs = find(subs);
            end
            if isa(x, 'ADI')
                x.jac{ix} = x.jac{ix}(:,subs);
            end
        end
        
        function x = replaceADI(x, y, replace, subs)
            if islogical(subs)
                subs = find(subs);
            end
            if islogical(replace)
                replace = find(replace);
            end
            if isa(x, 'ADI')
                x.val(replace) = y.val(subs);
                for i = 1 : numel(x.jac)
                    [n, m] = size(x.jac{i});
                    if n == m
                        x.jac{i}(replace,replace) = y.jac{i}(subs,subs);
                    elseif m == numel(subs)
                        x.jac{i}(replace,subs) = y.jac{i}(subs,subs);
                    end
                end
            else
                x(replace) = y(subs);
            end
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

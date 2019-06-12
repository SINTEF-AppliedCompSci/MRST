classdef PhaseMixingCoefficientsLV < StateFunction
    properties
        useCompactEvaluation = true;
    end
    
    methods
        function gp = PhaseMixingCoefficientsLV(model, varargin)
            gp@StateFunction(model, varargin{:});
        end

        function v = evaluateOnDomain(prop, model, state)
            eos = model.EOSModel;
            [p, T, x, y] = model.getProps(state,...
                'pressure', 'temperature', 'liquidMoleFractions', 'vaporMoleFractions');
            acf = eos.fluid.acentricFactors;
            
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
            
            L = struct('Si', {Si_L}, 'A', {A_L}, 'B', {B_L}, 'Bi', {Bi}, 'A_ij', {A_ij});
            V = struct('Si', {Si_V}, 'A', {A_V}, 'B', {B_V}, 'Bi', {Bi}, 'A_ij', {A_ij});
            if model.water
                v = {struct(); L; V};
            else
                v = {L; V};
            end
        end
    end
end
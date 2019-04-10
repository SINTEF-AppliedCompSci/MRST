classdef PhaseMixingCoefficientsLV < GridProperty
    properties
    end
    
    methods
        function gp = PhaseMixingCoefficientsLV(model, varargin)
            gp@GridProperty(model, varargin{:});
        end

        function v = evaluateOnDomain(prop, model, state)
            eos = model.EOSModel;
            [p, T, x, y] = model.getProps(state,...
                'pressure', 'temperature', 'liquidMoleFractions', 'vaporMoleFractions');
            acf = eos.fluid.acentricFactors;
            
            [A_ij, Bi] = eos.getMixingParameters(p, T, acf, iscell(x));
            [Si_L, A_L, B_L] = eos.getPhaseMixCoefficients(x, A_ij, Bi);
            [Si_V, A_V, B_V] = eos.getPhaseMixCoefficients(y, A_ij, Bi);
            
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
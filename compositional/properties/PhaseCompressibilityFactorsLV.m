classdef PhaseCompressibilityFactorsLV < StateFunction
    properties
    end
    
    methods
        function gp = PhaseCompressibilityFactorsLV(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PhaseMixingCoefficients', 'ComponentPhaseMoleFractions'});
        end

        function v = evaluateOnDomain(prop, model, state)
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
            
            Z_L = eos.computeCompressibilityZ(p, x, L_mix.A, L_mix.B, L_mix.Si, L_mix.Bi, true);
            Z_V = eos.computeCompressibilityZ(p, y, V_mix.A, V_mix.B, V_mix.Si, V_mix.Bi, false);
            
            s = getSampleAD(p, x{:}, y{:});
            
            Z_L = model.AutoDiffBackend.convertToAD(Z_L, s);
            Z_V = model.AutoDiffBackend.convertToAD(Z_V, s);

            if isfield(state, 'cellJacMap')
                arg = {state.cellJacMap};
            else
                arg = {};
            end
            Z_L = eos.setZDerivatives(Z_L, L_mix.A, L_mix.B, arg{:});
            Z_V = eos.setZDerivatives(Z_V, V_mix.A, V_mix.B, arg{:});
            
            v = cell(1, nph);
            v{L_ix} = Z_L;
            v{V_ix} = Z_V;
            if model.water
                v{1} = ones(numelValue(Z_L));
            end
        end
    end
end
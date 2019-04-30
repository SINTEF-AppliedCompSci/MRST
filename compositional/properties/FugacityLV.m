classdef FugacityLV < GridProperty
    properties
    end
    
    methods
        function gp = FugacityLV(model, varargin)
            gp@GridProperty(model, varargin{:});
            gp = gp.dependsOn({'PhaseMixingCoefficients', 'ComponentPhaseMoleFractions', 'PhaseCompressibilityFactors'});
        end

        function f = evaluateOnDomain(prop, model, state)
            eos = model.EOSModel;
            p = model.getProps(state, 'pressure');
            [mix, mf, Z] = prop.getEvaluatedDependencies(state, 'PhaseMixingCoefficients', 'ComponentPhaseMoleFractions', 'PhaseCompressibilityFactors');
            
            ncomp = numel(eos.fluid.names);
            f = cell(ncomp, 2);
            wat = model.water;
            for i = 1:2
                xy = mf((1+wat):end, i + wat)';
                m = mix{i+wat};
                
                f(:, i) = model.EOSModel.computeFugacity(p, xy, Z{i+wat}, m.A, m.B, m.Si, m.Bi);
            end
        end
    end
end
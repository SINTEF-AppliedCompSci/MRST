classdef ComponentTotalFluxMolecularDiffusion < StateFunction
    % Total component flux including advection + molecular diffusion:
    %   v_i = sum_alpha v_{i,alpha} + J_i^{mol,diff}

    methods
        function sf = ComponentTotalFluxMolecularDiffusion(model, varargin)
            sf@StateFunction(model, varargin{:});
            sf = sf.dependsOn('ComponentPhaseFlux');
            sf = sf.dependsOn('ComponentTotalMolecularDiffFlux');
            sf.label = 'V_i';
        end

        function v = evaluateOnDomain(sf, model, state)
            ncomp = model.getNumberOfComponents();
            nph   = model.getNumberOfPhases();

            phaseFlux = sf.getEvaluatedDependencies(state, 'ComponentPhaseFlux');
            Jdiff     = sf.getEvaluatedDependencies(state, 'ComponentTotalMolecularDiffFlux');

            v = cell(ncomp, 1);
            for c = 1:ncomp
                vc = 0;
                for ph = 1:nph
                    f = phaseFlux{c, ph};
                    if ~isempty(f)
                        vc = vc + f;
                    end
                end
                if ~isempty(Jdiff) && numel(Jdiff) >= c && ~isempty(Jdiff{c})
                    vc = vc + Jdiff{c};
                end
                v{c} = vc;
            end
        end
    end
end
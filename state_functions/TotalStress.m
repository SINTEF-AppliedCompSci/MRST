classdef TotalStress < StateFunction
    % (Nonlinear)absolute permeability arising from poroelastic effects
    properties
    end
    
    methods
        function gp = TotalStress(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'pressure', 'pressure_matrix'}, 'state');
            gp = gp.dependsOn({'EffectiveStress'});
        end
        function v = evaluateOnDomain(prop, model, state)
            [p, pm] = model.getProps(state, 'pressure', 'pressure_matrix');
            effective_stress = prop.getEvaluatedDependencies(state, 'EffectiveStress');
            cM = model.mechModel.constitutive_coefficients_object;
            if model.G.griddim == 2
                pI = bsxfun(@times, p, [1, 1, 0]);
                pmI = bsxfun(@times, pm, [1, 1, 0]);
            else
                pI = bsxfun(@times, p, [1, 1, 1, 0, 0, 0]);
                pmI = bsxfun(@times, pm, [1, 1, 1, 0, 0, 0]);
            end
            v = effective_stress - cM.B_m.*pmI - cM.B_f.*pI;  
        end
    end
end
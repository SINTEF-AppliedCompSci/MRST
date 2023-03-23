classdef PoroelasticMatrixPoro < StateFunction
    % Effective poroelastic matrix porosity 
    properties
    end
    
    methods
        function gp = PoroelasticMatrixPoro(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'pressure', 'pressure_matrix', 'xd'}, 'state');
        end
        function m_poro = evaluateOnDomain(prop, model, state)
            % Get effective matrix porosity given changes in strain and 
            % matrix and fracture pressures. We assume a linear form where
            % refrence pressures and strain are zero. 
            cM = model.mechModel.constitutive_coefficients_object;
            ref_poro_m = model.rock_matrix.poro;
            [p, pm, xd] = model.getProps(state, 'pressure', 'pressure_matrix',...
                                                'xd');
            mechTerm = model.computeStrainTerms(xd);
            m_poro = ref_poro_m + (mechTerm.matrix + cM.invN_m.*pm...
                                   + cM.invQ.*p); % mechTerm.matrix = b_m:E
        end
    end
end
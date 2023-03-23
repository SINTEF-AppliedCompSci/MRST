classdef PoroelasticFracturePoro < StateFunction
    % Effective poroelastic fracture porosity 
    properties
    end
    
    methods
        function gp = PoroelasticFracturePoro(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'pressure', 'pressure_matrix', 'xd'}, 'state');
        end
        function poro_f = evaluateOnDomain(prop, model, state)
            % Get effective fracture porosity given changes in strain and 
            % matrix and fracture pressures. We assume a linear form where
            % refrence pressures and strain are zero.
            cM = model.mechModel.constitutive_coefficients_object;
            ref_poro_f = model.rock.poro;
            [p, pm, xd] = model.getProps(state, 'pressure', 'pressure_matrix',...
                                                'xd');
            mechTerm = model.computeStrainTerms(xd);
            poro_f = ref_poro_f + (mechTerm.fracture + cM.invQ.*pm...
                                   + cM.invN_f.*p); % mechTerm.fracture = b_f:E
        end
    end
end
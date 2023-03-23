classdef PoroelasticFracturePerm < StateFunction
    % (Nonlinear)absolute permeability arising from poroelastic effects
    properties
    end
    
    methods
        function gp = PoroelasticFracturePerm(model, varargin)
            gp@StateFunction(model, varargin{:});
            if isfield(model.rock, 'nonLinearPerm')
                gp = gp.dependsOn({'PoroelasticFracturePoro'});
            end
        end
        function perm = evaluateOnDomain(prop, model, state)
            r = model.rock;
            if isfield(r, 'nonLinearPerm')
                poro_f = prop.getEvaluatedDependencies(state, 'PoroelasticFracturePoro');
                perm = prop.evaluateFunctionOnDomainWithArguments(r.nonLinearPerm, poro_f);
            end
        end
    end
end
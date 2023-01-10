classdef WellborePhaseUpwindFlag < StateFunction
    
    methods
        %-----------------------------------------------------------------%
        function cf = WellborePhaseUpwindFlag(model)

            cf@StateFunction(model);
%             cf = cf.dependsOn('ComponentPhaseFlux');

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function flag = evaluateOnDomain(prop, model, state)
            
%             v = prop.getEvaluatedDependencies(state,...
%                         'ComponentPhaseFlux');
            % @@TODO: Implement explicit Brennier & Jaffré upwinding
            vt   = value(state.massFlux);
            nc   = model.getNumberOfComponents();
            nph  = model.getNumberOfPhases();
            flag = repmat({vt >= 0}, 1, nph);
%             flag = cellfun(@(v) value(v) >= 0, v, 'UniformOutput', false);
            
        end
        %-----------------------------------------------------------------%
        
    end
    
end
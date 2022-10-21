classdef WellboreSurfaceRate < StateFunction
    
    methods
        %-----------------------------------------------------------------%
        function cf = WellboreSurfaceRate(model)

            cf@StateFunction(model);
            cf = cf.dependsOn({'ComponentPhaseFlux', 'PhaseUpwdinFlag'}, 'FlowPropertyFunctions');
            cf = cf.dependsOn({'SurfaceDensity'}, 'PVTpropertyFunctions');
            cf.label = 'q_s';

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function qs = evaluateOnDomain(prop, model, state)
            
            [q, rhoS, flag] = model.getProps(state, ...
                'ComponentPhaseFlux', 'SurfaceDensity', 'PhaseUpwindFlag');
            [~, iif] = model.getInletSegments();
            rhoS = cellfun(@(flag, rho) ...
                model.parentModel.operators.faceUpstr(flag, rho), ...
                flag, rhoS, 'UniformOutput', false);
            qs = cellfun(@(q, rho) q(iif)./rho(iif), ...
                    q, rhoS, 'UniformOutput', false);
            qs = cellfun(@(qs) model.parentModel.operators.fluxSum(qs), ...
                    qs, 'UniformOutput', false);
            
        end
        %-----------------------------------------------------------------%
        
    end
    
end
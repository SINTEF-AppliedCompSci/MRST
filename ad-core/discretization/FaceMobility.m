classdef FaceMobility < GridProperty & UpwindProperty
    properties
        
    end
    
    methods
        function fm = FaceMobility(backend, upwinding)
            fm@GridProperty(backend);
            fm@UpwindProperty(upwinding)
            fm = fm.dependsOn('PhaseUpwindFlag');
            fm = fm.dependsOn('Mobility', 'FlowPropertyFunctions');
        end
        
        function fmob = evaluateOnDomain(prop, model, state)
            flag = prop.getEvaluatedDependencies(state, 'PhaseUpwindFlag');
            mob = model.getProp(state, 'Mobility');
            nph = numel(mob);
            fmob = cell(1, nph);
            for i = 1:nph
                fmob{i} = prop.faceUpstream(state, flag{i}, mob{i});
            end
        end
    end
end
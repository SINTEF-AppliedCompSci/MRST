classdef FaceComponentMobility < GridProperty & UpwindProperty
    properties
        
    end
    
    methods
        function fm = FaceComponentMobility(backend, upwinding)
            fm@GridProperty(backend);
            fm@UpwindProperty(upwinding)
            fm = fm.dependsOn('ComponentMobility', 'FlowPropertyFunctions');
            fm = fm.dependsOn('PhaseUpwindFlag');
        end
        
        function mobf = evaluateOnDomain(prop, model, state)
            flag = prop.getEvaluatedDependencies(state, 'PhaseUpwindFlag');
            mob = model.getProps(state, 'ComponentMobility');
            [ncomp, nph] = size(mob);
            mobf = cell(ncomp, nph);
            for c = 1:ncomp
                for ph = 1:nph
                    if ~isempty(mob{c, ph})
                        mobf{c, ph} = prop.faceUpstream(state, flag{ph}, mob{c, ph});
                    end
                end
            end
        end
    end
end
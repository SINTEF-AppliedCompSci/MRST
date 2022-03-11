classdef FaceComponentMobilityAverage < StateFunction
    properties
        
    end
    
    methods
        function fm = FaceComponentMobilityAverage(backend, upwinding)
            fm@StateFunction(backend);
            fm = fm.dependsOn('ComponentMobility', 'FlowPropertyFunctions');
        end
        
        function mobf = evaluateOnDomain(prop, model, state)
            mob = model.getProps(state, 'ComponentMobility');
            [ncomp, nph] = size(mob);
            mobf = cell(ncomp, nph);
            for c = 1:ncomp
                for ph = 1:nph
                    if ~isempty(mob{c, ph})
                        mobf{c, ph} = model.operators.faceAvg(mob{c, ph});
                    end
                end
            end
        end
    end
end
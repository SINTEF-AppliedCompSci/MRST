classdef FaceMobilityAverage < StateFunction
    properties
        
    end
    
    methods
        function fm = FaceMobilityAverage(backend, upwinding)
            fm@StateFunction(backend);
            fm = fm.dependsOn('Mobility', 'FlowPropertyFunctions');
        end
        
        function fmob = evaluateOnDomain(prop, model, state)
            mob = model.getProp(state, 'Mobility');
            nph = numel(mob);
            fmob = cell(1, nph);
            for i = 1:nph
                fmob{i} = model.operators.faceAvg(mob{i});
            end
        end
    end
end
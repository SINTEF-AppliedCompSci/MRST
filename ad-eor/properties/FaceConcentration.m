classdef FaceConcentration < StateFunction & UpwindProperty
    properties
        
    end
    
    methods
        function gp = FaceConcentration(backend, upwinding)
            gp@StateFunction(backend);
            gp@UpwindProperty(upwinding)
            gp = gp.dependsOn('PhaseUpwindFlag');
        end
        
        function fcp = evaluateOnDomain(prop, model, state)
            flag = prop.getEvaluatedDependencies(state, 'PhaseUpwindFlag');
            cp = model.getProps(state, 'polymer');
            fcp = prop.faceUpstream(model, state, flag{1}, cp);            
        end
    end
end
classdef FaceMobility < GridProperty & UpwindProperty
    properties
        
    end
    
    methods
        function fm = FaceMobility(backend, upwinding)
            fm@GridProperty(backend);
            fm@UpwindProperty(upwinding)
        end
        
        function fmob = evaluateOnDomain(prop, model, state)
            [mob, flag] = model.getProps(state, 'Mobility', 'PhaseUpwindFlag');
            nph = numel(mob);
            fmob = cell(1, nph);
            for i = 1:nph
                fmob{i} = prop.faceUpstream(state, flag{i}, mob{i});
            end
        end
    end
end
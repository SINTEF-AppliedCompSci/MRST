classdef FaceMobility < GridProperty & UpwindProperty
    properties
        
    end
    
    methods
        function fm = FaceMobility(backend, upwinding)
            fm@GridProperty(backend);
            fm@UpwindProperty(upwinding)
        end
        
        function fmob = evaluateOnDomain(prop, model, state)
            [mob, pot] = model.getProps(state, 'Mobility', 'PhasePotentialDifference');
            nph = numel(mob);
            fmob = cell(1, nph);
            for i = 1:nph
                flag = pot{i} <= 0;
                fmob{i} = prop.faceUpstream(state, flag, mob{i});
            end
        end
    end
end
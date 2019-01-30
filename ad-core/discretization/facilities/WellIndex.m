classdef WellIndex < GridProperty
    properties

    end
    
    methods

        function WI = evaluateOnDomain(prop, model, state)
            map = model.getProp(state, 'FacilityWellMapping');
            wellSol = state.wellSol;
            W = model.getWellStruct(map.active);
            cstatus = vertcat(wellSol(map.active).cstatus);
            WI = vertcat(W.WI).*cstatus;
        end
    end
end
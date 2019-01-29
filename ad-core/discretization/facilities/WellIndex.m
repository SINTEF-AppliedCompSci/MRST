classdef WellIndex < GridProperty
    properties

    end
    
    methods

        function WI = evaluateOnDomain(prop, model, state)
            wellSol = state.wellSol;
            active = model.getIndicesOfActiveWells(wellSol);
            W = model.getWellStruct(active);
            cstatus = vertcat(wellSol(active).cstatus);
            WI = vertcat(W.WI).*cstatus;
        end
    end
end
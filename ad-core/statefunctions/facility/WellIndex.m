classdef WellIndex < StateFunction
    % Get the well index (productivity or injectivity index, equivialent to
    % a transmissibility between the reservoir and the wellbore).
    properties

    end
    
    methods
        function gp = WellIndex(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn('FacilityWellMapping');
        end
        function WI = evaluateOnDomain(prop, model, state)
            map = prop.getEvaluatedDependencies(state, 'FacilityWellMapping');
            wellSol = state.wellSol;
            W = model.getWellStruct(map.active);
            cstatus = vertcat(wellSol(map.active).cstatus);
            WI = vertcat(W.WI).*cstatus;
        end
    end
end
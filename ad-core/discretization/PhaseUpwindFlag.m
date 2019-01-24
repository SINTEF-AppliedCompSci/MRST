classdef PhaseUpwindFlag < GridProperty
    properties

    end
    
    methods
        function flags = evaluateOnDomain(prop, model, state)
            pot = model.getProp(state, 'PhasePotentialDifference');
            nph = numel(pot);
            flags = cell(1, nph);
            for i = 1:nph
                flags{i} = pot{i} <= 0;
            end
        end
    end
end
classdef PhaseUpwindFlag < GridProperty
    properties

    end
    
    methods
        function gp = PhaseUpwindFlag(varargin)
            gp@GridProperty(varargin{:});
            gp = gp.dependsOn('PhasePotentialDifference');
        end

        function flags = evaluateOnDomain(prop, model, state)
            pot = prop.getEvaluatedDependencies(state, 'PhasePotentialDifference');
            nph = numel(pot);
            flags = cell(1, nph);
            for i = 1:nph
                flags{i} = pot{i} <= 0;
            end
        end
    end
end
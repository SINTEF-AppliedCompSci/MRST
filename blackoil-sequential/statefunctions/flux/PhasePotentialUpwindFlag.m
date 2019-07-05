classdef PhasePotentialUpwindFlag < StateFunction
    properties

    end
    
    methods
        function gp = PhasePotentialUpwindFlag(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn('flux', 'state');
            gp = gp.dependsOn('Mobility', 'FlowPropertyFunctions');
            gp = gp.dependsOn({'PhaseInterfacePressureDifferences', 'Transmissibility', 'TotalFlux'});
        end

        function flags = evaluateOnDomain(prop, model, state)
            [G, T, vT] = prop.getEvaluatedDependencies(state, ...
                            'PhaseInterfacePressureDifferences', 'Transmissibility', 'TotalFlux');
            mob = model.getProp(state, 'Mobility');
            flag_array = multiphaseUpwindIndices(G, vT, T, mob, model.operators.faceUpstr);
            nph = numel(mob);
            flags = cell(1, nph);
            for i = 1:nph
                flags{i} = flag_array(:, i);
            end
        end
    end
end
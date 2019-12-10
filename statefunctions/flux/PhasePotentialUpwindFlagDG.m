classdef PhasePotentialUpwindFlagDG < PhasePotentialUpwindFlag
    properties
    end
    
    methods
        function gp = PhasePotentialUpwindFlagDG(varargin)
            gp@PhasePotentialUpwindFlag(varargin{:});
        end

        function flags = evaluateOnDomain(prop, model, state)
            mob = model.getProp(state, 'Mobility');
            nph = numel(mob);
            flags = cell(1, nph);
            if strcmpi(state.type, 'face')
                [G, T, vT] = prop.getEvaluatedDependencies(state, ...
                                'PhaseInterfacePressureDifferences', 'Transmissibility', 'TotalFlux');
                if prop.includeTotalVelocity
                    v = vT;
                else
                    v = zeros(size(value(vT)));
                end
                flag_array = multiphaseUpwindIndices(G, v, T, mob, model.operators.faceUpstr);
                for i = 1:nph
                    flags{i} = flag_array(:, i);
                end
            else
                [flags{:}] = deal(true);
            end
        end
    end
end
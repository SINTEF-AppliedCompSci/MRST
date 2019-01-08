classdef PhaseFlux < GridProperty
    properties

    end
    
    methods
        function v = evaluateOnGrid(prop, model, state)
            [mob, kgrad] = model.getProps(state, 'FaceMobility', 'PermeabilityPotentialGradient');
            nph = numel(mob);
            v = cell(1, nph);
            for i = 1:nph
                v{i} = -mob{i}.*kgrad{i};
            end
        end
    end
end
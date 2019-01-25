classdef PhaseFlux < GridProperty & UpwindProperty
    properties

    end
    
    methods
        function fm = PhaseFlux(backend, upwinding)
            fm@GridProperty(backend);
            fm@UpwindProperty(upwinding)
        end

        
        function v = evaluateOnDomain(prop, model, state)
            [mob, kgrad, flag] = model.getProps(state, ...
                'Mobility', 'PermeabilityPotentialGradient', 'PhaseUpwindFlag');
            nph = numel(mob);
            v = cell(1, nph);
            for i = 1:nph
                mobf = prop.faceUpstream(state, flag{i}, mob{i});
                v{i} = -mobf.*kgrad{i};
            end
        end
    end
end
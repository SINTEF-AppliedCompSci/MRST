classdef WellPhaseFlux < GridProperty
    properties

    end
    
    methods

        function q_ph = evaluateOnDomain(prop, model, state)
            map = model.getProp(state, 'FacilityWellMapping');
            W = map.W;
            [dp, wi] = model.getProps(state, 'PerforationPressureGradient', 'WellIndex');
            mob = model.ReservoirModel.getProps(state, 'Mobility');
            
            mobw = cellfun(@(x) x(map.cells), mob, 'UniformOutput', false);
            nph = numel(mob);
            
            isInj = dp < 0;
            if any(isInj)
                compi = vertcat(W.compi);
                compi_perf = compi(map.perf2well, :);
                mobt = zeros(sum(isInj), 1);
                for i = 1:nph
                    mobt = mobt + mobw{i}(isInj);
                end
                for i = 1:nph
                    mobw{i}(isInj) = mobt.*compi_perf(isInj, i);
                end
            end
            q_ph = cell(1, nph);
            for i = 1:nph
                q_ph{i} = mobw{i}.*wi.*dp;
            end
        end
    end
end
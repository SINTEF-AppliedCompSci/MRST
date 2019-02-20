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
            
            isInjector = map.isInjector(map.perf2well);
            
            Tdp = -wi.*dp;
            injection = Tdp > 0;
            crossflow = injection & ~isInjector;
            if any(injection)
                compi = vertcat(W.compi);
                if any(crossflow)
                    % warning('Crossflow occuring in %d perforations', sum(crossflow));
                    % Compute cross flow for this phase. The approach here
                    % is to calculate (as doubles) the volumetric inflow of
                    % all phases into the well-bore. If a well has
                    % cross-flow, the phase distribution of the
                    % cross-flowing volume is assumed to reflect the inflow
                    % conditions, neglecting density change throughout the
                    % wellbore.
                    q_wb = bsxfun(@times, value(mobw), value(Tdp));
                    q_wb = -min(q_wb, 0);
                    nw = numel(W);
                    comp = zeros(nw, nph);
                    for i = 1:nph
                        comp(:, i) = accumarray(map.perf2well, q_wb(:, i));
                    end
                    % Normalize to get volume fractions
                    compT = sum(comp, 2);
                    comp = bsxfun(@rdivide, comp, compT);
                    % We computed this everywhere, we only need it where
                    % crossflow is occuring.
                    keep = compT > 0 & ~map.isInjector;
                    compi(keep, :) = comp(keep, :);
                end
                compi_perf = compi(map.perf2well, :);
                mobt = zeros(sum(injection), 1);
                for i = 1:nph
                    mobt = mobt + mobw{i}(injection);
                end
                for i = 1:nph
                    mobw{i}(injection) = mobt.*compi_perf(injection, i);
                end
            end
            q_ph = cell(1, nph);
            for i = 1:nph
                q_ph{i} = mobw{i}.*Tdp;
            end
%             
%             if any(injection)
%                 compi
%                 value(q_ph)
%             end
        end
    end
end
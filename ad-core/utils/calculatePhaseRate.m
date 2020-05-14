function q_ph = calculatePhaseRate(Tdp, mobw, map, allowCrossFlow, model)
    isInjector = map.isInjector(map.perf2well);
    W = map.W;
    vTdp = value(Tdp);
    injection = vTdp > 0;
    production = ~injection & vTdp ~= 0;
    crossflow = (injection & ~isInjector) | ...
                (production & isInjector);
    crossflow = crossflow & allowCrossFlow;
    
    nph = numel(mobw);
    if any(injection)
        compi = vertcat(W.compi);
        if any(crossflow)
            dispif(model.verbose > 1, 'Crossflow occuring in %d perforations\n', sum(crossflow));
            % Compute cross flow for this phase. The approach here
            % is to calculate (as doubles) the volumetric inflow of
            % all phases into the well-bore. If a well has
            % cross-flow, the phase distribution of the
            % cross-flowing volume is assumed to reflect the inflow
            % conditions, neglecting density change throughout the
            % wellbore.
            q_wb = bsxfun(@times, value(mobw), vTdp);
            compi = crossFlowMixture(q_wb, compi, map);
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
end
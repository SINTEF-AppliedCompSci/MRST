function tot_strap = computeTotalStructuralTrapping(Gt, ta, rock, seainfo)

    surface_pressure = 0;
    press_deviation = 0;

    
    % Get trapping capacity:
    [~, strap, ~, ~, ~] = compute_trapcap(Gt, ta, rock, seainfo, ...
        surface_pressure, 'press_deviation',press_deviation);

    
    % total structural trapping only:
    tot_strap = sum(strap)/giga; % Mt

    
    % an assertion:
    ttrap_cap = [];
    num_traps = numel(unique(ta.trap_regions))-1;
    for tix=1:num_traps
        ttrap_cap(tix) = sum(strap(ta.traps == tix)); % kg
    end
    trap_cap = ttrap_cap/giga; % Mt
    assert(sum(trap_cap) - tot_strap < sqrt(eps));
    
end
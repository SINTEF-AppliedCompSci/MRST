function Ma = forecastCurve(initState, states, model, ta)

    % Using each instance of "states" (including initial) to make forecast curve:
    
    % NB: but if initial state didn't contain any CO2, don't forecast that
    % state, rather add forecast of 0 explicitly.
    if sum(initState.s(:,2)) == 0
        sts_base = {states{:}};
    else
        sts_base = {initState, states{:}};
    end
    Ma = zeros(numel(sts_base),1);
    for i = 1:numel(sts_base)
        Ma(i) = massAtInfinity( model.G, model.rock, sts_base{i}.pressure, ...
            sts_base{i}.s(:,2), sts_base{i}.sGmax, sts_base{i}.s(:,1), 0, ...
            model.fluid, ta, [], 'p_future',initState.pressure ) / 1e12; % Gt
    end
    if sum(initState.s(:,2)) == 0
        Ma = [0; Ma]; % now it includes forecast of initial state
    end

end
        
        
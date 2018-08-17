function initState = makeMyInitialState(Gt, seainfo)

    % Initial state is CO2-free.
    sG = zeros(Gt.cells.num,1);
    sF = 1 - sG;
    
    initState.s         = [sF, sG];
    initState.sGmax     = initState.s(:,2);
    
    gravity on;
    surface_pressure = 0;
    p_init = seainfo.water_density * norm(gravity()) * Gt.cells.z + surface_pressure;
    press_deviation = 0;
    p_init = p_init .* (1+press_deviation/100);
    initState.pressure  = p_init;

end
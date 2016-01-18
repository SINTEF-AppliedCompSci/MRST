function exploreOptWellNCS_postProcess( Gt, init, optim, other )
% Post-processing for optimizeFormations run using exploreOptWellNCS


    % Reconstruct the fluid structure which wasn't saved to avoid large
    % .mat files (~400 MB)
    other.fluid = reconstructFluidStructure( Gt, other.rock, other.opt, other.dh );

    reports_init  = makeReports_extras(Gt, {other.initState, init.states{:}}, ...
                                other.rock, other.fluid, init.schedule, ...
                                other.residual, other.traps, other.dh, ...
                                init.wellSols);


    reports_optim = makeReports_extras(Gt, {other.initState, optim.states{:}}, ...
                            other.rock, other.fluid, optim.schedule, ...
                            other.residual, other.traps, other.dh, ...
                            optim.wellSols);


    selectedResultsMultiplot(Gt, reports_optim, [2], ...
                             'plot_plume', false, ...
                             'plot_well_numbering', true, ...
                             'plot_distrib', false);


    %compareWellrates(init.schedule, optim.schedule, other.fluid.rhoGS);
    compareWellrates_viaWellSols(init.wellSols, optim.wellSols, ...
                                 init.schedule, other.fluid.rhoGS);


    % initial rates
    selectedResultsMultiplot(Gt, reports_init, [2 4 6], ...
                             'background', 'totalCO2', ...
                             'plot_traps', true);

    % optimized rates
    selectedResultsMultiplot(Gt, reports_optim, [1], ...
                             'background', 'totalCO2', ...
                             'plot_traps', true);


    % overpressure
    selectedResultsMultiplot(Gt, reports_optim, [61], ...
                             'background', 'overpressure', ...
                             'init_state', other.initState, ...
                             'plot_traps', true, 'plot_distrib', false, ...
                             'plot_well_numbering', true);

                     
    % Also:
    optim.states = computeOverpressure(other.initState, optim.states, optim.schedule);
    % figure; plotToolbar(Gt, {other.initState, optim.states{:}})
    % plotGrid(Gt, 'EdgeAlpha',0.1, 'FaceColor','none')
    % plotWellSols(optim.wellSols); set(gcf,'name','optim: bhp-wells')
    % 
    % 
    init.states = computeOverpressure(other.initState, init.states, init.schedule);
    % figure; plotToolbar(Gt, {other.initState, init.states{:}})
    % plotGrid(Gt, 'EdgeAlpha',0.1, 'FaceColor','none')
    % plotWellSols(init.wellSols); set(gcf,'name','init: bhp-wells')
    %     

    % sat plots just before migration period
    figure
    subplot(1,2,1)
    plotCellData(Gt, init.states{60}.s(:,2), 'EdgeAlpha',0.1)
    plotGrid(Gt, 'EdgeAlpha',0.1, 'FaceColor','none')
    colorbar
    axis equal tight
    title('Initial')

    subplot(1,2,2)
    plotCellData(Gt, optim.states{60}.s(:,2), 'EdgeAlpha',0.1)
    plotGrid(Gt, 'EdgeAlpha',0.1, 'FaceColor','none')
    colorbar
    axis equal tight
    title('Optimized')

    % press plots just before migration period (likely when the highest
    % overpressure occurred)
    figure
    subplot(1,2,1)
    plotCellData(Gt, init.states{60}.pressDev, 'EdgeAlpha',0.1)
    plotGrid(Gt, 'EdgeAlpha',0.1, 'FaceColor','none')
    colorbar
    axis equal tight
    title('Initial')

    subplot(1,2,2)
    plotCellData(Gt, optim.states{60}.pressDev, 'EdgeAlpha',0.1)
    plotGrid(Gt, 'EdgeAlpha',0.1, 'FaceColor','none')
    colorbar
    axis equal tight
    title('Optimized')

    seafloor_depth = other.opt.ref_depth; % based on naming convention in optimizeFormation
    water_density = other.opt.rhoW; % based on naming convention in optimizeFormation
    [~, P_limit] = computeOverburdenPressure(Gt, other.rock, ...
        seafloor_depth, ...
        water_density);
    plotFormationPressureChanges(optim.states, other.initState.pressure, ...
        'P_lim', P_limit)



end


function exploreOptWellNCS_postProcess( Gt, init, optim, other, varargin )
% Post-processing for optimizeFormations run using exploreOptWellNCS

% NB: ensure experimental/project/.../SampleProps2D and /CO2props is not on
% path, otherwise errors will occur. Use corresponding files in co2lab/...

    moduleCheck('mrst-gui')
    moduleCheck('ad-core')

    opt.savePlots   = false;
    opt.fmName      = [];
    opt.figDirName  = [];
    
    opt.plotWells           = false;
    opt.plotWellRates       = true;
    opt.plotInventory       = false;
    opt.plotPressureChanges = false;
    opt.plotOther           = false;
    
    opt.plotSatFields       = false;
    opt.SatFieldStates      = []; % can be an array
    
    opt.plotPressFields     = false;
    opt.PressFieldStates    = []; % can be an array
    
    opt.warningOn = true;
    
    opt = merge_options(opt, varargin{:});

    
    %% Reconstruct the fluid structure which wasn't saved to avoid large
    % .mat files (~400 MB)
    if ~isfield(other,'fluid')
        other.fluid = reconstructFluidStructure( Gt, other.rock, other.opt, other.dh );
    end
    
    %% Make Initial and Optimized Reports
    reports_init  = makeReports_extras(Gt, {other.initState, init.states{:}}, ...
                                other.rock, other.fluid, init.schedule, ...
                                other.residual, other.traps, other.dh, ...
                                init.wellSols, 'warningOn',opt.warningOn);


    reports_optim = makeReports_extras(Gt, {other.initState, optim.states{:}}, ...
                            other.rock, other.fluid, optim.schedule, ...
                            other.residual, other.traps, other.dh, ...
                            optim.wellSols, 'warningOn',opt.warningOn);

    %% Plot wells in formation
    if opt.plotWells
        selectedResultsMultiplot(Gt, reports_optim, [2], ...
            'plot_plume', false, 'plot_well_numbering', true, ...
            'plot_distrib', false);
    end

    %% Plot well rates, compute total injected/leaked, print to table  
    if opt.plotWellRates 
        %compareWellrates(init.schedule, optim.schedule, other.fluid.rhoGS);
        num_optim_wells = compareWellrates_viaWellSols(init.wellSols, optim.wellSols, ...
                                     init.schedule, other.fluid.rhoGS, Gt, other.traps);

        [total_inj, total_leaked] = totals_Injected_and_Leaked(Gt, other.rock, ...
                       other.fluid, optim.states, optim.schedule, optim.wellSols, ...
                       false);
                   
        seainfo = getSeaInfo(other.opt.modelname, other.opt.refRhoCO2); % contains seainfo.dis_max
        fmCap   = getTrappingInfo(Gt, other.rock, seainfo, 'mapPlotOn',false);
        
        % full storage potentials:
        strapCap = fmCap.breakdown.structural_trapping_capacity; % Gt
        rtrapCap = fmCap.breakdown.residual_trapping_capacity;   % Gt
        dtrapCap = fmCap.breakdown.dissolved_trapping_capacity;  % Gt
        
        % achieved (simulated) storage (in kg):
        strap_achieved = reports_optim(end).masses(2) + reports_optim(end).masses(5) + reports_optim(end).masses(6); % resStruct + freeStruct + subTrap, kg
        rtrap_achieved = reports_optim(end).masses(3) + reports_optim(end).masses(4); % resTrap + freeRes, kg
        dtrap_achieved = reports_optim(end).masses(1); % resDis
        other_achieved = reports_optim(end).masses(7); % freeMov
        

%         fprintf('Formation       | StrapCap(Gt) | RtrapCap(Gt)| NumWellsPlaced | NumOptimalWells | Total inj. (Gt) | Percent leaked | Percent StrapAchieved | Percent RtrapAchieved | MovePlume(Gt) \n');
        fprintf('%16s&    %8.2f  &    %8.2f  &     %8.0f  &       %8.0f  &       %8.2f  &        %8.1f  &          %8.1f  &        %8.1f  &      %8.3f    \\\\ \n', ...
            other.opt.modelname, ...
            strapCap, ...
            rtrapCap, ...
            numel([init.schedule.control(1).W.cells]), ...
            num_optim_wells, ...
            convertTo(convertFrom(total_inj, mega), giga), ...
            (total_leaked/total_inj)*100, ...
            (convertTo(strap_achieved, giga*kilo)/strapCap)*100, ...
            (convertTo(rtrap_achieved, giga*kilo)/rtrapCap)*100, ...
            convertTo(other_achieved, giga*kilo) );
        
%         plotStorageBreakdownsPie( strapCap, rtrapCap, dtrapCap, ...
%             convertTo(strap_achieved,giga*kilo), convertTo(rtrap_achieved,giga*kilo), convertTo(dtrap_achieved,giga*kilo) )

        if opt.savePlots
           % 2, 3, 4, 5 figures generated. Saving only well rates.
           %saveas(gcf, [opt.figDirName '/' opt.fmName '_initVsOptimRates'], 'fig')
           drawnow
           pause(1)
           export_fig(figure(2), [opt.figDirName '/' opt.fmName '_initVsOptimRates'], '-png','-transparent')
        end
    end
    
    %% Plot trapping mechanisms inventory over time
    if opt.plotInventory
        
        % make inventory plots of init and optim:
        h = figure;

        % init inventory
        subplot(1,2,1); plot(1);
        ax = get(h, 'currentaxes');  
        % load all timesteps up to last plotted one (a bit of a hack)
        plotTrappingDistribution(ax, reports_init, 'legend_location', 'northeast');
        fsize = 24;
        set(get(gca, 'xlabel'), 'fontsize', fsize)
        set(get(gca, 'ylabel'), 'fontsize', fsize)
        set(gca,'fontsize', fsize);
        %set(gcf, 'position', [1, 1, 850, 850]);

        % optim inventory
        subplot(1,2,2); plot(1);
        ax = get(h, 'currentaxes');  
        % load all timesteps up to last plotted one (a bit of a hack)
        plotTrappingDistribution(ax, reports_optim, 'legend_location', 'northeast');
        fsize = 24;
        set(get(gca, 'xlabel'), 'fontsize', fsize)
        set(get(gca, 'ylabel'), 'fontsize', fsize)
        set(gca,'fontsize', fsize);
        set(gcf, 'position', [1, 1, 1835 788]);

        if opt.savePlots
            %saveas(gcf, [opt.figDirName '/' opt.fmName _inventory'], 'fig')
            drawnow
            pause(1)
            export_fig(gcf, [opt.figDirName '/' opt.fmName '_inventory'], '-png','-transparent')
        end
    end    
    
    %% Plot Saturation Fields and/or Pressure Fields
    if opt.plotSatFields
        
        % initial rates
        selectedResultsMultiplot(Gt, reports_init, opt.SatFieldStates, ...
                                 'background', 'totalCO2', ...
                                 'plot_traps', true);

        % optimized rates
        selectedResultsMultiplot(Gt, reports_optim, opt.SatFieldStates, ...
                                 'background', 'totalCO2', ...
                                 'plot_traps', true);
    end
    
    if opt.plotPressFields
        % overpressure
        selectedResultsMultiplot(Gt, reports_optim, opt.PressFieldStates, ...
                                 'background', 'overpressure', ...
                                 'init_state', other.initState, ...
                                 'plot_traps', true, 'plot_distrib', false, ...
                                 'plot_well_numbering', true);
    end
    
    %% Compute pressure field when highest overpressure occurred
    if opt.plotOther
        
        fprintf('Optimized Rates:\n')
        [optim.states, maxPressDev_stepNum] = computeOverpressure(other.initState, optim.states, optim.schedule);
        % figure; plotToolbar(Gt, {other.initState, optim.states{:}})
        % plotGrid(Gt, 'EdgeAlpha',0.1, 'FaceColor','none')
        % plotWellSols(optim.wellSols); set(gcf,'name','optim: bhp-wells')
        % 
        % 
        fprintf('Initial Rates:\n')
        [init.states, ~] = computeOverpressure(other.initState, init.states, init.schedule);
        % figure; plotToolbar(Gt, {other.initState, init.states{:}})
        % plotGrid(Gt, 'EdgeAlpha',0.1, 'FaceColor','none')
        % plotWellSols(init.wellSols); set(gcf,'name','init: bhp-wells')
        %     

        % sat plots at end of migration
        figure
        subplot(1,2,1)
        plotCellData(Gt, init.states{end}.s(:,2), 'EdgeAlpha',0.1)
        plotGrid(Gt, 'EdgeAlpha',0.1, 'FaceColor','none')
        colorbar
        axis equal tight
        title('Initial')

        subplot(1,2,2)
        plotCellData(Gt, optim.states{end}.s(:,2), 'EdgeAlpha',0.1)
        plotGrid(Gt, 'EdgeAlpha',0.1, 'FaceColor','none')
        colorbar
        axis equal tight
        title('Optimized')

        % press plots when the highest overpressure occurred
        figure
        subplot(1,2,1)
        plotCellData(Gt, init.states{maxPressDev_stepNum}.pressDev, 'EdgeAlpha',0.1)
        plotGrid(Gt, 'EdgeAlpha',0.1, 'FaceColor','none')
        colorbar
        axis equal tight
        title('Initial')

        subplot(1,2,2)
        plotCellData(Gt, optim.states{maxPressDev_stepNum}.pressDev, 'EdgeAlpha',0.1)
        plotGrid(Gt, 'EdgeAlpha',0.1, 'FaceColor','none')
        colorbar
        axis equal tight
        title('Optimized')
    
    end
    
    %% Plot Pressure Changes
    if opt.plotPressureChanges
        
        seafloor_depth = other.opt.ref_depth; % based on naming convention in optimizeFormation
        water_density = other.opt.rhoW; % based on naming convention in optimizeFormation
        [P_over, ~] = computeOverburdenPressure(Gt, other.rock, ...
            seafloor_depth, ...
            water_density);
        plotFormationPressureChanges(optim.states, other.initState.pressure, ...
            P_over, Gt, optim.schedule);
        % use res for a summary display
        
        if opt.savePlots
            % several plots generated. Saving last one.
            %saveas(gcf, [opt.figDirName '/' opt.fmName _inventory'], 'fig')
            drawnow
            pause(1)
            export_fig(gcf, [opt.figDirName '/' opt.fmName '_pressureChanges'], '-png','-transparent')
        end
    end

end

function [total_inj, total_leaked] = totals_Injected_and_Leaked(Gt, rock2D, ...
    fluid, states, schedule, wellSols, printOutput)
% Calculations are based on leak_penalizer function written in
% optimizeRates.m


%     % take only wellSols.qGr corresponding to injection period (control step 1)
%     wellSols = wellSols( schedule.step.control == 1 );
%     
%     % compute total co2 mass injected
%     rates           = zeros(numel(wellSols), nw);
%     vols            = zeros(numel(wellSols), nw);
%     timeStepSizes   = schedule.step.val(schedule.step.control == 1);
%     for i = 1:numel(wellSols)
%         rates(i,:) = [wellSols{i}.qGr]; % m3/s
%         vols(i,:) = [wellSols{i}.qGr] .* timeStepSizes(i); % m3
%     end
%     total_vols_inj  = sum(vols); % totals per well
%     totalVol        = sum(total_vols_inj); % total injected into formation
%     totalMass       = totalVol * co2RefRho; % m3 * kg/m3 = kg
    
    
    % compute total co2 mass leaked by end of simulation
    numSteps = numel(states);
    tSteps = (1:numSteps)';
    dts = schedule.step.val;
    
    vol_inj = 0; % @@
    for step = 1:numSteps
        sol = wellSols{tSteps(step)};
        state = states{tSteps(step)}; %@@ +1?
        nW = numel(sol);
        pBHP = zeros(nW, 1); % place holder
        qGs = vertcat(sol.qGs);
        qWs = vertcat(sol.qWs);
        p = state.pressure;
        %p = compute_hydrostatic_pressure(model.G, rho_water, surf_press);
        sG = state.s(:,2);
        dt = dts(step);
        injInx = (vertcat(sol.sign) > 0);
        
        % calculate "M^inj"
        vol_inj  = vol_inj +  dt * spones(ones(1, nW)) * (injInx .* qGs); % m3
        mass_inj = vol_inj * fluid.rhoGS / 1e9; % Mt
        
        if (tSteps(step) == numSteps)
            % calculate "M^accum"
            bG         = fluid.bG(p);
            pvol       = Gt.cells.volumes .* Gt.cells.H .* rock2D.poro;
            if isfield(rock2D,'ntg')
                pvol   = Gt.cells.volumes .* Gt.cells.H .* rock2D.poro .* rock2D.ntg;
                % accounting for possible net-to-gross data
            end
            vol_accum  = ones(1, Gt.cells.num) * (pvol .* fluid.pvMultR(p) .* bG .* sG); % m3
            mass_accum = vol_accum * fluid.rhoGS / 1e9; % Mt
            
            if printOutput
                fprintf('Total injected: %f (m3) or %f (Mt)\n', vol_inj, mass_inj);
                fprintf('Total leaked: %f (m3) or %f (Mt)\n', (vol_inj - vol_accum), (mass_inj - mass_accum)); 
            end
        end
    end
    
    % Pass out results:
    total_inj       = mass_inj;
    total_leaked    = mass_inj - mass_accum;
    
    


end


function exploreOptWellNCS_postProcess( Gt, init, optim, other, varargin )
% Post-processing for optimizeFormations run using exploreOptWellNCS

% NB: ensure experimental/project/.../SampleProps2D and /CO2props is not on
% path, otherwise errors will occur. Use corresponding files in co2lab/...

    %close all
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
    opt.outputOn = true;
    
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
    
    %% Plot Pressure Changes
    if opt.plotPressureChanges
        
        seafloor_depth = other.opt.ref_depth; % based on naming convention in optimizeFormation
        water_density = other.opt.rhoW; % based on naming convention in optimizeFormation
        [P_over, ~] = computeOverburdenPressure(Gt, other.rock, ...
            seafloor_depth, ...
            water_density);
        res = plotFormationPressureChanges(optim.states, other.initState.pressure, ...
            P_over, Gt, optim.schedule, 'outputOn', opt.outputOn);
        % use res for a summary display
        
        if opt.savePlots
            % several plots generated. Saving last one.
            %saveas(gcf, [opt.figDirName '/' opt.fmName _inventory'], 'fig')
            drawnow
            pause(1)
            export_fig(gcf, [opt.figDirName '/' opt.fmName '_pressureChanges'], '-png','-transparent')
        end
    end

    %% Plot well rates, compute total injected/leaked, print to table  
    if opt.plotWellRates 
        %compareWellrates(init.schedule, optim.schedule, other.fluid.rhoGS);
        num_optim_wells = compareWellrates_viaWellSols(init.wellSols, optim.wellSols, ...
                                     init.schedule, other.fluid.rhoGS, Gt, other.traps, ...
                                     init, optim);

        [total_inj, total_leaked] = totals_Injected_and_Leaked(Gt, other.rock, ...
                       other.fluid, optim.states, optim.schedule, optim.wellSols, ...
                       false);
                   
        seainfo = getSeaInfo(other.opt.modelname, other.opt.refRhoCO2); % contains seainfo.dis_max
        fmCap   = getTrappingInfo(Gt, other.rock, seainfo, 'mapPlotOn',false);
        
        % full storage potentials:
        strapCap = fmCap.breakdown.structural_trapping_capacity; % Gt
        rtrapCap = fmCap.breakdown.residual_trapping_capacity;   % Gt
        dtrapCap = fmCap.breakdown.dissolved_trapping_capacity;  % Gt
        totCap   = fmCap.breakdown.total_trapping_capacity;      % Gt
        
        pv = Gt.cells.volumes .* Gt.cells.H .* other.rock.poro; % pore volume (m3)
        if isfield(other.rock,'ntg')
            pv = pv .* other.rock.ntg;
        end
        
        % achieved (simulated) storage (in kg):
        strap_achieved = reports_optim(end).masses(2) + reports_optim(end).masses(5) + reports_optim(end).masses(6); % resStruct + freeStruct + subTrap, kg
        rtrap_achieved = reports_optim(end).masses(3) + reports_optim(end).masses(4); % resTrap + freeRes, kg
        dtrap_achieved = reports_optim(end).masses(1); % resDis
        other_achieved = reports_optim(end).masses(7); % freeMov
        
        balance_error = (convertTo((strap_achieved + rtrap_achieved + dtrap_achieved + other_achieved),giga*kilo) - (total_inj - total_leaked));
        assert( abs(balance_error) < 1e-12 )

% %         fprintf('Formation       | StrapCap(Gt) | RtrapCap(Gt) | DtrapCap(Gt) | Total inj. (Gt) | Leaked (Gt,Perc.) | Perc.FracPress | StrapAchieved (Gt,Perc.) | RtrapAchieved (Gt,Perc.) | DtrapAch (Gt,Perc.) | MovePlume(Gt) \n');
%         %fprintf('%16s&   %8.2f  &   %8.2f  &   %8.2f  &    %8.2f   & %8.2f (%3.1f) &    %3.1f    &  %8.2f (%3.1f)  &   %8.2f (%3.1f)   &  %8.2f (%3.1f)  &  %8.2f    \\\\ \n', ...
%         fprintf('%16s&   %8.3f  &   %8.3f  &   %8.3f  &    %8.3f   & %8.3f (%3.2f) &    %3.2f    &  %8.3f (%3.2f)  &   %8.3f (%3.2f)   &  %8.3f (%3.2f)  &  %8.3f    \\\\ \n', ...
%             other.opt.modelname, ...
%             strapCap, ...
%             rtrapCap, ... %numel([init.schedule.control(1).W.cells]), num_optim_wells, 
%             dtrapCap, ...
%             total_inj, ...
%             total_leaked, (total_leaked/total_inj)*100, ...
%             res.worst_percent_of_fracPress_reached, ...
%             convertTo(strap_achieved, giga*kilo), (convertTo(strap_achieved, giga*kilo)/strapCap)*100, ...
%             convertTo(rtrap_achieved, giga*kilo), (convertTo(rtrap_achieved, giga*kilo)/rtrapCap)*100, ...
%             convertTo(dtrap_achieved, giga*kilo), (convertTo(dtrap_achieved, giga*kilo)/dtrapCap)*100, ...    
%             convertTo(other_achieved, giga*kilo) );
        
        % compute (simulated) storage efficiency:
        rhoG_init   = other.fluid.rhoG(other.initState.pressure); % kg/m3
        full_fm_cap = sum(rhoG_init .* pv)/1e12; % Gt
        Seff        = (total_inj - total_leaked) / full_fm_cap;
        
        % simplier table output:
        if strcmpi(other.opt.btype,'pressure')
        %fprintf('Formation       | Total inj. (Gt) | Leaked (Gt,Perc.) | Cap (Gt) | Seff (Perc.) | Perc.FracPress | StrapAchieved (Gt,Perc.) | RtrapAchieved (Gt,Perc.) | MovePlume(Gt) \n');
        fprintf('%16s&      %8.3f   &   %8.3f (%3.2f) & %8.3f &        %3.2f  &       %3.2f &     %8.3f (%3.2f)  &        %8.3f (%3.2f)   &  %8.3f    \\\\ \n', ...
            other.opt.modelname, ...
            total_inj, ...
            total_leaked, (total_leaked/total_inj)*100, ...
            total_inj - total_leaked, ...
            Seff * 100, ...
            res.worst_percent_of_fracPress_reached, ...
            convertTo(strap_achieved, giga*kilo), (convertTo(strap_achieved, giga*kilo)/strapCap)*100, ...
            convertTo(rtrap_achieved, giga*kilo), (convertTo(rtrap_achieved, giga*kilo)/rtrapCap)*100, ...  
            convertTo(other_achieved, giga*kilo) );
        
        elseif strcmpi(other.opt.btype,'flux')
        %fprintf('Formation       | Total inj. (Gt) | Seff (Perc.) | Perc.FracPress | StrapAchieved (Gt,Perc.) | RtrapAchieved (Gt,Perc.) | MovePlume(Gt) \n');
        fprintf('%16s&      %8.3f   &         %3.2f  &       %3.2f &     %8.3f (%3.2f)  &        %8.3f (%3.2f)   &  %8.3f    \\\\ \n', ...
            other.opt.modelname, ...
            total_inj, ...
            Seff * 100, ...
            res.worst_percent_of_fracPress_reached, ...
            convertTo(strap_achieved, giga*kilo), (convertTo(strap_achieved, giga*kilo)/strapCap)*100, ...
            convertTo(rtrap_achieved, giga*kilo), (convertTo(rtrap_achieved, giga*kilo)/rtrapCap)*100, ...  
            convertTo(other_achieved, giga*kilo) );
            
        end
        
        
        % only trapping capacity output:
%         fprintf('Formation       | Vp (km3) | StrapCap(Gt) (Perc) | RtrapCap(Gt) (Perc) | DtrapCap(Gt) (Perc) | Total (Gt) \n');
%         fprintf('%16s&   %8.2f  & %8.3f  (%3.2f) &   %8.3f (%3.2f) &   %8.3f (%3.2f) &    %8.3f  \\\\ \n', ...
%             other.opt.modelname, ...
%             convertTo(sum(pv), kilo*kilo*kilo), ...
%             strapCap, (strapCap/totCap)*100, ...
%             rtrapCap, (rtrapCap/totCap)*100, ... %numel([init.schedule.control(1).W.cells]), num_optim_wells, 
%             dtrapCap, (dtrapCap/totCap)*100, ...
%             totCap);

%         plotStorageBreakdownsPie( strapCap, rtrapCap, dtrapCap, ...
%             convertTo(strap_achieved,giga*kilo), convertTo(rtrap_achieved,giga*kilo), convertTo(dtrap_achieved,giga*kilo) )

        if opt.savePlots
           % 2, 3, 4, 5 figures generated. Saving only well rates.
           %saveas(gcf, [opt.figDirName '/' opt.fmName '_initVsOptimRates'], 'fig')
           drawnow
           pause(1)
           export_fig(figure(3), [opt.figDirName '/' opt.fmName '_initVsOptimRates'], '-png','-transparent')
        end
    end
    
    %% Plot trapping mechanisms inventory over time
    if opt.plotInventory
        
        % make inventory plots of init and optim:
        h = figure; set(gcf,'Position',[3690 814 1428 545])

        % init inventory
        subplot(1,2,1); plot(1);
        ax = get(h, 'currentaxes');  
        % load all timesteps up to last plotted one (a bit of a hack)
        plotTrappingDistribution(ax, reports_init, 'legend_location', 'northeast');
        fsize = 20;
        set(get(gca, 'xlabel'), 'fontsize', fsize)
        set(get(gca, 'ylabel'), 'fontsize', fsize)
        set(gca,'fontsize', fsize);
        %set(gcf, 'position', [1, 1, 850, 850]);

        % optim inventory
        subplot(1,2,2); plot(1);
        ax = get(h, 'currentaxes');  
        % load all timesteps up to last plotted one (a bit of a hack)
        plotTrappingDistribution(ax, reports_optim, 'legend_location', 'northeast');
        fsize = 20;
        set(get(gca, 'xlabel'), 'fontsize', fsize)
        set(get(gca, 'ylabel'), 'fontsize', fsize)
        set(gca,'fontsize', fsize);
        %set(gcf, 'position', [1, 1, 1835 788]);

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
    

end

function [total_inj, total_leaked] = totals_Injected_and_Leaked(Gt, rock2D, ...
    fluid, states, schedule, wellSols, printOutput)
% Calculations are based on leak_penalizer function written in
% optimizeRates.m

% NEED to UPDATE if DISSOLUTION WAS SIMULATED


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
        mass_inj = vol_inj * fluid.rhoGS / 1e12; % Gt
        
        if (tSteps(step) == numSteps)
            % calculate "M^accum"
            bG         = fluid.bG(p);
            pvol       = Gt.cells.volumes .* Gt.cells.H .* rock2D.poro;
            if isfield(rock2D,'ntg')
                pvol   = Gt.cells.volumes .* Gt.cells.H .* rock2D.poro .* rock2D.ntg;
                % accounting for possible net-to-gross data
            end
            vol_accum  = ones(1, Gt.cells.num) * (pvol .* fluid.pvMultR(p) .* bG .* sG); % m3
            mass_accum = vol_accum * fluid.rhoGS / 1e12; % Gt
            
            if printOutput
                fprintf('Total injected: %f (m3) or %f (Gt)\n', vol_inj, mass_inj);
                fprintf('Total leaked: %f (m3) or %f (Gt)\n', (vol_inj - vol_accum), (mass_inj - mass_accum)); 
            end
        end
    end
    
    % Pass out results:
    total_inj       = mass_inj;                 % Gt
    total_leaked    = mass_inj - mass_accum;    % Gt
    
    


end


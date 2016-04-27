function Seff = exploreOptWellNCS_postProcess( Gt, init, optim, other, varargin )
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
    opt.plotWellsPlaced     = true;
    opt.plotWellsRemaining  = true;
    opt.plotInventory       = true;
    opt.plotPressureChanges = true;
    opt.plotOther           = false;
    opt.plot_co2_height     = true;
    opt.plot_fraction_overburden = true;
    
    opt.arbitrary_well_cost = []; % USD/well
    
    opt.warningOn = true;
    opt.outputOn = true;
    
    opt = merge_options(opt, varargin{:});
    
    % check that p_lim_factor is a variable, otherwise it was likely not
    % used during simulation. In that case, set it here to avoid break
    if ~isfield(other.opt,'p_lim_factor')
       other.opt.p_lim_factor = 0.9; 
    end

    
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
                        
    %% Plot saturation (or co2 height) field
    if opt.plot_co2_height
        last_inj_step = numel(optim.schedule.step.val(optim.schedule.step.control == 1)) + 1;
        last_mig_step = numel(reports_init);
        if ~isfield(other.rock,'ntg')
           other.rock.ntg = 1; 
        end
        selectedResultsMultiplot(Gt, reports_init, [last_inj_step last_mig_step], ...
                             'background', 'totalCO2', ...
                             'background_threshold', 0.01, ...
                             'ref_temp', other.opt.ref_temp, ...
                             'ref_depth', other.opt.ref_depth, ...
                             'temp_grad', other.opt.temp_grad, ...
                             'poro', other.rock.poro, ...
                             'ntg', other.rock.ntg, ...
                             'plot_traps', true, ...
                             'plot_well_numbering', true, ...
                             'plume_threshold', 0.3, ...
                             'plot_distrib', false);
        % adjust axis of each subfigure in current figure
        hfig = gcf; set(hfig,'Position',[3765 178 655 547]);
        axs = findobj(hfig.Children,'type','axes');
        for i=1:numel(axs)
           set(hfig,'CurrentAxes',axs(i)); axis equal tight
        end

        selectedResultsMultiplot(Gt, reports_optim, [last_inj_step last_mig_step], ...
                             'background', 'totalCO2', ...
                             'background_threshold', 0.01, ...
                             'ref_temp', other.opt.ref_temp, ...
                             'ref_depth', other.opt.ref_depth, ...
                             'temp_grad', other.opt.temp_grad, ...
                             'poro', other.rock.poro, ... 
                             'ntg', other.rock.ntg, ...
                             'plot_traps', true, ...
                             'plot_well_numbering', true, ...
                             'plume_threshold', 0.3, ...
                             'plot_distrib', false);
        % adjust axis of each subfigure in current figure
        hfig = gcf; set(hfig,'Position',[4432 179 655 547]);
        axs = findobj(hfig.Children,'type','axes');
        for i=1:numel(axs)
           set(hfig,'CurrentAxes',axs(i)); axis equal tight
        end
    end
                     
    %% Plot Pmax/Pover at a given time step
    if opt.plot_fraction_overburden
        [P_over, ~] = computeOverburdenPressure(Gt, other.rock, other.opt.ref_depth, other.opt.rhoW);
        selectedResultsMultiplot(Gt, reports_optim, [last_inj_step last_mig_step], ...
                             'background', 'fraction_overburden_reached', ...
                             'overburden', P_over, ...
                             'plot_plume', false, ...
                             'background_threshold', [], ...
                             'plot_traps', true, ...
                             'plot_well_numbering', true, ...
                             'plot_distrib', false);
        % adjust axis of each subfigure in current figure
        hfig = gcf; set(hfig,'Position',[4432 2 655 546]);
        axs = findobj(hfig.Children,'type','axes');
        for i=1:numel(axs)
           set(hfig,'CurrentAxes',axs(i)); axis equal tight
        end
    end

    %% Plot wells in formation
    if opt.plotWells
        selectedResultsMultiplot(Gt, reports_optim, [2], ...
            'plot_plume', false, 'plot_well_numbering', true, ...
            'plot_distrib', false);
    end
    
    %% Plot Pressure Changes
    [P_over, ~] = computeOverburdenPressure(Gt, other.rock, other.opt.ref_depth, other.opt.rhoW);
    res = plotFormationPressureChanges(optim.states, other.initState.pressure, ...
        P_over, P_over*other.opt.p_lim_factor, ...
        Gt, optim.schedule, 'outputOn', opt.outputOn, ...
        'makePlot', opt.plotPressureChanges);
    % use res for a summary display
        
    if opt.plotPressureChanges && opt.savePlots
        % several plots generated. Saving last one.
        %saveas(gcf, [opt.figDirName '/' opt.fmName _inventory'], 'fig')
        drawnow
        pause(1)
        export_fig(gcf, [opt.figDirName '/' opt.fmName '_pressureChanges'], '-png','-transparent')
    end
        

    %% Plot well rates, compute total injected/leaked, print to table  
    % extend for when optimal rates obtained by iteratively increasing cp
    if opt.plotWellRates 
        %compareWellrates(init.schedule, optim.schedule, other.fluid.rhoGS);
        if ~other.opt.account_well_cost
        num_optim_wells = compareWellrates_viaWellSols(init.wellSols, optim.wellSols, ...
                                     init.schedule, other.fluid.rhoGS, Gt, other.traps, ...
                                     init, optim, other, 'plotWellsPlaced', opt.plotWellsPlaced, ...
                                     'plotWellsRemaining', opt.plotWellsRemaining);
        else
        num_optim_wells = compareWellrates_viaWellSols(init.wellSols, optim.wellSols, ...
                                     init.schedule, other.fluid.rhoGS, Gt, other.traps, ...
                                     init, optim, other, 'plotWellsPlaced', opt.plotWellsPlaced, ...
                                     'plotWellsRemaining', opt.plotWellsRemaining, ...
                                     'plotCriticalRate', true, ...
                                     'arbitrary_well_cost', opt.arbitrary_well_cost);
        end
        if opt.savePlots
           % 2, 3, 4, 5 figures generated. Saving only well rates.
           %saveas(gcf, [opt.figDirName '/' opt.fmName '_initVsOptimRates'], 'fig')
           drawnow
           pause(1)
           export_fig(figure(3), [opt.figDirName '/' opt.fmName '_initVsOptimRates'], '-png','-transparent')
        end
    end

    % this is the total leaked by simulation
    [total_inj, total_leaked] = totals_Injected_and_Leaked(Gt, other.rock, ...
                   other.fluid, optim.states, optim.schedule, optim.wellSols, ...
                   false);
               
    % this is the total leaked by prediction
    model.G = Gt;
    model.rock = other.rock;
    model.fluid = other.fluid;
    [~, Mi_tot, Ma] = leak_penalizer_at_infinity_Rerun(model, optim.wellSols, optim.states, optim.schedule, ...
        other.opt.leak_penalty, other.opt.surface_pressure, other.opt.rhoW, other.traps, 'plotsOn',true);
    total_leaked_future = Mi_tot(end) - Ma(end); % Gt
    
    % to compare simulated leakage to predicted leakage
    % injected mass | remaining | leaked (%) | predicted remaining | predicted leaked (%)
    fprintf('%5.4f & %5.4f & %5.4f (%3.3f) & %5.4f & %5.4f (%3.3f) \\\\ \n', ...
        total_inj*1e3, ...
        (total_inj - total_leaked)*1e3, ...
        total_leaked*1e3, (total_leaked/total_inj)*100, ...
        (total_inj - total_leaked_future)*1e3, ...
        total_leaked_future*1e3, (total_leaked_future/total_inj)*100 );
    

    % Get maximal retaining capacity (trapping breakdown). Then,
    % determine how much of each trapping has been achieved. (Use final
    % pressure state to determine maximal retaining capacity, or at the
    % very least, surface pressure should be included in calculation of
    % hydrostatic pressure).
    if ~strcmpi(other.opt.modelname,'Synthetic')
        seainfo = getSeaInfo(other.opt.modelname, other.opt.refRhoCO2); % contains seainfo.dis_max
    else
        seainfo = getSeaInfo('NorthSea', other.opt.refRhoCO2);
    end
    %fmCap   = getTrappingInfo(Gt, other.rock, seainfo, 'mapPlotOn',false, ...
    %                           'surf_press', other.opt.surface_pressure);
    fmCap   = getTrappingInfo(Gt, other.rock, seainfo, 'mapPlotOn',false, ...
                                'press_field', optim.states{end}.pressure); % @@ compare impact on open bdrys

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
    %masses = max([resDis, resStruc, resTrap, freeRes, freeStruc, subtrap, freeMov], 0);
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
%             res.worst_percent_of_OBpress_reached, ...
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
    %fprintf('Formation       | Total inj. (Gt) | Leaked (Gt,Perc.) | Cap (Gt) | Seff (Perc.) | Perc.FracPress | StrapAchieved (Gt,Perc.) | RtrapAchieved (Gt,Perc.) | MovePlume(Gt) | cp \n');
    fprintf('%16s&      %8.3f   &   %8.3f (%3.2f) & %8.3f &        %3.2f  &       %3.2f &     %8.3f (%3.2f)  &        %8.3f (%3.2f)   &  %8.3f  &  %1.0e   \\\\ \n', ...
        other.opt.modelname, ...
        total_inj, ...
        total_leaked, (total_leaked/total_inj)*100, ...
        total_inj - total_leaked, ...
        Seff * 100, ...
        res.worst_percent_of_OBpress_reached, ...
        convertTo(strap_achieved, giga*kilo), (convertTo(strap_achieved, giga*kilo)/strapCap)*100, ...
        convertTo(rtrap_achieved, giga*kilo), (convertTo(rtrap_achieved, giga*kilo)/rtrapCap)*100, ...  
        convertTo(other_achieved, giga*kilo), ...
        other.opt.pressure_penalty );
    
    % Accounting for well cost:
    % NB: if qGs >= 2*sqrt(eps), well treated as active
    % # wells placed | total_inj(Mt) | total_inj - total_leaked |
    % total_leaked (total_leaked/total_inj)*100 | S(MillionUSD) | # active wells | 
    % S - # active wells * well_initial_cost
    if ~isempty(opt.arbitrary_well_cost)
    num_active_wells = numel([optim.schedule.control(1).W.val] >= 2*sqrt(eps)); % @@ test threshold
    alpha = 1;
    if isfield(other.opt,'alpha')
        alpha = other.opt.alpha;
    end
    fprintf('%2.0f & %5.4f & %5.4f & %5.4f (%3.3f) & %6.3f & %2.0f & %2.2f & %6.3f \\\\ \n', ...
        numel([optim.schedule.control(1).W.val]), ...
        total_inj*1e3, ...
        (total_inj - total_leaked)*1e3, ...
        total_leaked*1e3, (total_leaked/total_inj)*100, ...
        optim.obj_val_total/1e6, ... % when obj val is in USD
        num_active_wells, ...
        alpha, ...
        (optim.obj_val_total - num_active_wells*opt.arbitrary_well_cost)/1e6);
    % generates figure 12
    account_well_cost_rerun(Gt, other.fluid, optim.wellSols, optim.schedule, ...
        other.opt.well_initial_cost, other.opt.well_operation_cost, ...
        other.opt.co2_tax_credit, alpha, other.opt.nonlinear_well_cost);
    end

    elseif strcmpi(other.opt.btype,'flux')
    %fprintf('Formation       | Total inj. (Gt) | Seff (Perc.) | Perc.FracPress | StrapAchieved (Gt,Perc.) | RtrapAchieved (Gt,Perc.) | MovePlume(Gt) | cp \n');
    fprintf('%16s&      %8.3f   &         %3.2f  &       %3.2f &     %8.3f (%3.2f)  &        %8.3f (%3.2f)   &  %8.3f    &  %1.0e \\\\ \n', ...
        other.opt.modelname, ...
        total_inj, ...
        Seff * 100, ...
        res.worst_percent_of_OBpress_reached, ...
        convertTo(strap_achieved, giga*kilo), (convertTo(strap_achieved, giga*kilo)/strapCap)*100, ...
        convertTo(rtrap_achieved, giga*kilo), (convertTo(rtrap_achieved, giga*kilo)/rtrapCap)*100, ...  
        convertTo(other_achieved, giga*kilo), ...
        other.opt.pressure_penalty );

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


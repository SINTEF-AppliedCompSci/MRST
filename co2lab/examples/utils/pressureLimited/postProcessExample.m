function postProcessExample(Gt, init, optim, other)
% Post-processing of Bjarmeland pressure-limited example    

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    if ~isfield(other.rock,'ntg')
       other.rock.ntg = 1; 
    end

    
    %% Well placement and trapping structure
    figure; clf;
    wcells = [init.schedule.control(1).W.cells];
    mapPlot(gcf, Gt, 'traps', other.ta.traps, ...
        'trapcolor', [0.5 0.5 0.5], 'trapalpha', 0.7, ...
        'rivers', other.ta.cell_lines, 'rivercolor', [0 0 1], ...
        'maplines', 40, 'wellcells', wcells, 'well_numbering', true);
    colorizeCatchmentRegions(Gt, other.ta);
    plotGrid(Gt, 'FaceColor','none','EdgeAlpha',0.1)
    axis equal tight off
    title('Well placement & trapping structure')
    
    
    %% Make Initial and Optimized Reports
    reports_init  = makeReports(Gt, {other.initState, init.states{:}}, ...
                                    other.rock, other.fluid, init.schedule, ...
                                    other.residual, other.ta, []);


    reports_optim = makeReports(Gt, {other.initState, optim.states{:}}, ...
                                other.rock, other.fluid, optim.schedule, ...
                                other.residual, other.ta, []);

    %% Plot co2 saturation
    % NB: 'totalCO2' will plot tons of mass of CO2 per lateral square meter
    last_inj_step = numel(optim.schedule.step.val(optim.schedule.step.control == 1)) + 1;
    last_mig_step = numel(reports_init);

    selectedResultsMultiplot(Gt, reports_init, [last_inj_step last_mig_step], ...
                         'background', 'totalCO2', ...
                         'background_threshold', 0.01, ...
                         'backgroundalpha', 1, ...
                         'ref_temp', other.ref_temp, ...
                         'ref_depth', other.ref_depth, ...
                         'temp_grad', other.temp_grad, ...
                         'poro', other.rock.poro, ...
                         'ntg', other.rock.ntg, ...
                         'plot_traps', true, ...
                         'plot_well_numbering', true, ...
                         'plume_threshold', 0.3, ...
                         'plot_distrib', false, ...
                         'trapmethod', [], ...
                         'ta', other.ta);

    % adjust axis of each subfigure in current figure
    hfig = gcf; set(hfig,'Position',[1 1 1000 700]);
    axs = findobj(hfig.Children,'type','axes');
    for i=1:numel(axs)
       set(hfig,'CurrentAxes',axs(i)); axis equal tight
       plotFaces(Gt, boundaryFaces(Gt));
    end


    selectedResultsMultiplot(Gt, reports_optim, [last_inj_step last_mig_step], ...
                         'background', 'totalCO2', ...
                         'background_threshold', 0.01, ...
                         'backgroundalpha', 1, ...
                         'ref_temp', other.ref_temp, ...
                         'ref_depth', other.ref_depth, ...
                         'temp_grad', other.temp_grad, ...
                         'poro', other.rock.poro, ...
                         'ntg', other.rock.ntg, ...
                         'plot_traps', true, ...
                         'plot_well_numbering', true, ...
                         'plume_threshold', 0.3, ...
                         'plot_distrib', false, ...
                         'trapmethod', [], ...
                         'ta', other.ta);

    % adjustments
    hfig = gcf; set(hfig,'Position',[1 1 1000 700]);
    axs = findobj(hfig.Children,'type','axes');
    for i=1:numel(axs)
       set(hfig,'CurrentAxes',axs(i)); axis equal tight
       plotFaces(Gt, boundaryFaces(Gt));
    end
   

    %% Max fraction Overburden pressure reached

    % Overburden pressure (recomputed to be sure)
    P_over = computeOverburdenPressure(Gt, other.rock, other.ref_depth, other.rhoW);
    selectedResultsMultiplot(Gt, reports_optim, [2], ...
                         'background', P_over, ...
                         'plot_plume', false, ...
                         'background_threshold', 0, ...
                         'backgroundalpha', 0.8, ...
                         'plot_traps', false, ...
                         'plot_well_numbering', true, ...
                         'plot_distrib', false, ...
                         'trapmethod', [], ...
                         'ta', other.ta );
    % adjustments
    title('Overburden pressure', 'fontsize',12)
    hfig = gcf; set(hfig,'Position',[1 1 500 700]);
    axs = findobj(hfig.Children,'type','axes');
    for i=1:numel(axs)
       set(hfig,'CurrentAxes',axs(i)); axis equal tight
    end
    
    
    % Max frac. reached                 
    tmp = [];
    for mx = 1:numel(reports_optim)
        tmp = [tmp, reports_optim(mx).sol.pressure];
    end
    maxP_encountered = max(tmp')';
    field = maxP_encountered./P_over;
                     
    selectedResultsMultiplot(Gt, reports_optim, [2], ...
                         'background', field, ...
                         'plot_plume', false, ...
                         'background_threshold', 0, ...
                         'backgroundalpha', 0.8, ...
                         'plot_traps', false, ...
                         'plot_well_numbering', true, ...
                         'plot_distrib', false, ...
                         'trapmethod', [], ...
                         'ta', other.ta );
    % adjustments
    hfig = gcf; set(hfig,'Position',[1 1 500 700]);
    axs = findobj(hfig.Children,'type','axes');
    for i=1:numel(axs)
       set(hfig,'CurrentAxes',axs(i)); axis equal tight
    end
    title({'Fraction overburden';'pressure reached'}); %, 'fontsize',12)
    

    %% Plot well rates
    compareWellrates(init.schedule, optim.schedule, other.fluid.rhoGS);
    legend('initial','optimized')


    %% Trapping inventory
    % includes forecast curve
    model.G = Gt;
    model.rock = other.rock;
    model.fluid = other.fluid;

    % forecast curves of initial and optim solution
    [~, ~, Ma_init] = leak_penalizer_at_infinity_Rerun(model, init.wellSols, ...
        init.states, init.schedule, other.leak_penalty, other.surface_pressure, ...
        other.rhoW, other.ta, 'plotsOn',true,'h',101);

    [~, ~, Ma_optim] = leak_penalizer_at_infinity_Rerun(model, optim.wellSols, ...
            optim.states, optim.schedule, other.leak_penalty, other.surface_pressure, ...
            other.rhoW, other.ta, 'plotsOn',true,'h',102);


    % make inventory plots of init and optim:
    h = figure; set(gcf,'Position',[1 1 1400 550])

    % init inventory
    subplot(1,2,1); plot(1);
    ax = get(h, 'currentaxes');  
    % load all timesteps up to last plotted one (a bit of a hack)
    plotTrappingDistribution(ax, reports_init, 'legend_location', 'northeast');
    fsize = 12;
    set(get(gca, 'xlabel'), 'fontsize', fsize)
    set(get(gca, 'ylabel'), 'fontsize', fsize)
    set(gca,'fontsize', fsize);
    % add forecast curve (Ma is in Gt)
    hold on;
    plot([0; cumsum(convertTo(init.schedule.step.val,year))], [0; Ma_init*1e3], ':b','LineWidth',5)
    legend off
    title('Initial')

    % optim inventory
    subplot(1,2,2); plot(1);
    ax = get(h, 'currentaxes');  
    % load all timesteps up to last plotted one (a bit of a hack)
    plotTrappingDistribution(ax, reports_optim, 'legend_location', 'northeast');
    fsize = 12;
    set(get(gca, 'xlabel'), 'fontsize', fsize)
    set(get(gca, 'ylabel'), 'fontsize', fsize)
    set(gca,'fontsize', fsize);
    % add forecast curve (Ma is in Gt)
    hold on;
    plot([0; cumsum(convertTo(optim.schedule.step.val,year))], [0; Ma_optim*1e3], ':b','LineWidth',5)

    % adjust yaxis of initial inventory to match yaxis of optim inventory
    hax = subplot(1,2,2);
    subplot(1,2,1);
    ylim([hax.YLim]);

    % add forecast curve label to legend
    hfig = gcf;
    str = findobj(hfig.Children,'type','legend');
    str = str.String;
    str(end+1) = {'Forecast'};
    subplot(1,2,2), hold on, legend(str);
    title('Optimized')


end

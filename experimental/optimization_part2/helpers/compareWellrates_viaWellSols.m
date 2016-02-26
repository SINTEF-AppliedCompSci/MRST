function num_optim_wells = compareWellrates_viaWellSols(wellSols1, wellSols2, schedule, co2RefRho, ...
    Gt, ta, init, optim)
%
% Updated to compute rates using qGs, not qGr.

    % take only wellSols.qGs corresponding to injection period (control step 1)
    
%     wS1 = cell2mat( wellSols1( schedule.step.control == 1 ) );
%     wS2 = cell2mat( wellSols2( schedule.step.control == 1 ) );
%     timeSizes = schedule.step.val(schedule.step.control == 1);
% 
%     % compute avg qGs per well (required if bhp-controlled, as qGs will
%     % vary over the injection period)
%     [~,nw] = size(wS1);
%     assert( numel(schedule.control(1).W) == nw )
%     for i=1:nw
%         
%         qGs1(i) = mean([wS1(:,i).qGs]); % m3/s
%         qGs2(i) = mean([wS2(:,i).qGs]);
%         % or keep track of max and min, and include as spread in bar-plot
%         
%         masses1(i) = mean([wS1(:,i).qGs]' .* co2RefRho .* timeSizes); % kg
%         masses2(i) = mean([wS2(:,i).qGs]' .* co2RefRho .* timeSizes); % kg
%     end
% 
%     %
%     avg_init_rates = qGs1; % m3/s
%     avg_opt_rates  = qGs2;
%     avg_init_masses = masses1; % kg
%     avg_optim_masses = masses2;
    
    %% 
    assert( strcmpi(init.schedule.control(1).W(1).type, 'rate') )
    init_rates = [init.schedule.control(1).W.val]; % should be the same as qGs if rate-wells
    opt_rates = [optim.schedule.control(1).W.val];
    
    injSec = cumsum(optim.schedule.step.val(optim.schedule.step.control == 1));
    injSec = injSec(end);
    tot_inj_init = sum([init.schedule.control(1).W.val] .* co2RefRho .* injSec ./ 1e9); % Mt
    tot_inj_opt = sum([optim.schedule.control(1).W.val] .* co2RefRho .* injSec ./ 1e9); % Mt
    
    %% Plot average rates
    val = [init_rates; opt_rates]';
    
    % values are currently in m3/sec using the CO2 reference density.  We
    % convert it to Mt/year
    val = val * co2RefRho * year / 1e6 / 1e3;
    
    % if only 1 well exists, perform hack to make bar plot work properly:
    if numel(init_rates) == 1
        val = [val; 0 0];
    end
    
    % determine how many figures are currently open, and avoid overwritting
    hfig = gcf;
    x = hfig.Number;
    figure(1+x); set(gcf,'Position',[2562 938 560 420])
    clf
    bar(val);
    set(gca, 'xlim', [0 numel(init_rates)+1])
    if numel(init_rates) == 1
        set(gca, 'xlim', [0 numel(init_rates)+1], 'XTick', [1:numel(init_rates) []]);
    end
    xlabel('Well', 'fontsize', 16);
    ylabel('Rate (Mt/year)', 'fontsize', 16);
    set(gca, 'fontsize', 16); % increase fontsize on axes.
    box off;
    legend(['Initial (total injected: ',num2str(tot_inj_init), ' Mt)'],...
           ['Optimal (total injected: ',num2str(tot_inj_opt),' Mt)'])
    
    % an optimal well is considered as one which has a rate > 0.01 Mt/yr
    %num_optim_wells = numel(find(qGs2 > 0.01));
    % an optimal well is considered as one which has rate above the minimum
    % rate constraint. @@ what to do for bhp-controls?
    min_rate = sqrt(eps);
    num_optim_wells = numel(opt_rates(opt_rates > 2*min_rate)); % @@ or set threshold for neglegible rate
    
%     %% Plot relative contributions: use this to describe the relative
%     % importance of the wells.
%     relative_contribution_init = avg_init_masses./sum(avg_init_masses);
%     relative_contribution_optim = avg_optim_masses./sum(avg_optim_masses);
%     figure(2+x);
%     clf
%     bar([relative_contribution_init; relative_contribution_optim]'.*100)
%     set(gca, 'xlim', [0 nw+1])
%     if nw == 1
%         set(gca, 'xlim', [0 nw+1], 'XTick', [1:nw []]);
%     end
%     xlabel('Well', 'fontsize', 16);
%     ylabel('Contribution to total injected mass (%)', 'fontsize', 16);
%     set(gca, 'fontsize', 16); % increase fontsize on axes.
%     box off;
%     legend(['Total injected mass: ',num2str(sum(avg_init_masses)/1e9), ' Mt'],...
%            ['Total injected mass: ',num2str(sum(avg_optim_masses)/1e9),' Mt'])
%     close
    
%     % Relative importance of wells
%     num_init_wells = nw;
%     num_optim_wells = numel(avg_optim_masses(avg_optim_masses > 2*min(avg_optim_masses)));
%     figure;
%     bar([relative_contribution_init./num_init_wells; relative_contribution_optim./num_optim_wells]')
%     xlabel('Well', 'fontsize', 16);
%     ylabel('Contribution to total injected mass (%)', 'fontsize', 16);
%     set(gca, 'fontsize', 16); % increase fontsize on axes.
%     box off;

    %% Plot mapPlot of wells and wells remaining:
    cinx_inj = [schedule.control(1).W.cells];
    figure(3+x); set(gcf,'Position',[2562 436 560 420])
    clf
    mapPlot(gcf, Gt, 'traps', ta.traps, ...
        'trapcolor', [0.5 0.5 0.5], 'trapalpha', 0.7, ...
        'rivers', ta.cell_lines, 'rivercolor', [1 0 0], ...
        'maplines', 20, 'wellcells', cinx_inj, 'well_numbering', true);
    colorizeCatchmentRegions(Gt, ta);
    plotGrid(Gt, 'FaceColor','none','EdgeAlpha',0.1)
    axis equal tight off
    title('Wells placed')
    
    
    figure(4+x); set(gcf,'Position',[3125 436 560 420])
    clf
    mapPlot(gcf, Gt, 'traps', ta.traps, ...
        'trapcolor', [0.5 0.5 0.5], 'trapalpha', 0.7, ...
        'rivers', ta.cell_lines, 'rivercolor', [1 0 0], ...
        'maplines', 20, 'wellcells', cinx_inj, 'well_numbering', true);
    colorizeCatchmentRegions(Gt, ta);
    plotGrid(Gt, 'FaceColor','none','EdgeAlpha',0.1)
    axis equal tight off
    
    % remove well number and dot of those wells deemed to have close to
    % zero-rate
    %masses_tmp = avg_optim_masses;
    %masses_tmp(avg_optim_masses < 2*min(avg_optim_masses)) = 0;
    rates_tmp = opt_rates;
    rates_tmp(opt_rates < 2*min_rate) = 0;
    hfig = gcf;
    wellLabels = findobj(hfig.Children.Children,'Type','Text');
    lines = findobj(hfig.Children.Children,'Type','Line');
    assert(numel(lines(end).XData) == numel(init_rates));
    Xdots_updated = lines(end).XData;
    Ydots_updated = lines(end).YData;
    rates_tmp_flipped = fliplr(rates_tmp); % ordering of wellLabels is backwards
    for i = 1:numel(rates_tmp) 
        if rates_tmp_flipped(i) == 0
            wellLabels(i).String='  '; 
        end
    end
    for i = 1:numel(rates_tmp) % ordering of well dots is forwards
        if rates_tmp(i) == 0
            Xdots_updated(i) = 0;
            Ydots_updated(i) = 0;
        end
    end
    lines(end).XData = Xdots_updated(Xdots_updated > 0);
    lines(end).YData = Ydots_updated(Ydots_updated > 0);
    title('Wells remaining')
    
    

end
function num_optim_wells = compareWellrates_viaWellSols(wellSols1, wellSols2, schedule, co2RefRho, ...
    Gt, ta)
%
% Updated to compute rates using qGs, not qGr.

    % take only wellSols.qGs corresponding to injection period (control step 1)
    
    wS1 = cell2mat( wellSols1( schedule.step.control == 1 ) );
    wS2 = cell2mat( wellSols2( schedule.step.control == 1 ) );
    timeSizes = schedule.step.val(schedule.step.control == 1);

    % compute avg qGs per well (required if bhp-controlled, as qGs will
    % vary over the injection period)
    [~,nw] = size(wS1);
    assert( numel(schedule.control(1).W) == nw )
    for i=1:nw
        
        qGs1(i) = mean([wS1(:,i).qGs]); % m3/s
        qGs2(i) = mean([wS2(:,i).qGs]);
        % or keep track of max and min, and include as spread in bar-plot
        
        masses1(i) = mean([wS1(:,i).qGs]' .* co2RefRho .* timeSizes); % kg
        masses2(i) = mean([wS2(:,i).qGs]' .* co2RefRho .* timeSizes); % kg
    end

    %
    avg_init_rates = qGs1; % m3/s
    avg_opt_rates  = qGs2;
    avg_init_masses = masses1; % kg
    avg_optim_masses = masses2;

    
    %% Plot average rates
    val = [avg_init_rates; avg_opt_rates]';
    
    % values are currently in m3/sec using the CO2 reference density.  We
    % convert it to Mt/year
    val = val * co2RefRho * year / 1e6 / 1e3;
    
    % if only 1 well exists, perform hack to make bar plot work properly:
    if nw == 1
        val = [val; 0 0];
    end
    
    figure; bar(val);
    set(gca, 'xlim', [0 nw+1])
    if nw == 1
        set(gca, 'xlim', [0 nw+1], 'XTick', [1:nw []]);
    end
    xlabel('Well', 'fontsize', 16);
    ylabel('Rate (Mt/year)', 'fontsize', 16);
    set(gca, 'fontsize', 16); % increase fontsize on axes.
    box off;
    
    % an optimal well is considered as one which has a rate > 0.01 Mt/yr
    %num_optim_wells = numel(find(qGs2 > 0.01));
    num_optim_wells = numel(avg_optim_masses(avg_optim_masses > 2*min(avg_optim_masses)));
    
    %% Plot relative contributions: use this to describe the relative
    % importance of the wells.
    relative_contribution_init = avg_init_masses./sum(avg_init_masses);
    relative_contribution_optim = avg_optim_masses./sum(avg_optim_masses);
    figure;
    bar([relative_contribution_init; relative_contribution_optim]'.*100)
    set(gca, 'xlim', [0 nw+1])
    if nw == 1
        set(gca, 'xlim', [0 nw+1], 'XTick', [1:nw []]);
    end
    xlabel('Well', 'fontsize', 16);
    ylabel('Contribution to total injected mass (%)', 'fontsize', 16);
    set(gca, 'fontsize', 16); % increase fontsize on axes.
    box off;
    legend(['Total injected mass: ',num2str(sum(avg_init_masses)/1e9), ' Mt'],...
           ['Total injected mass: ',num2str(sum(avg_optim_masses)/1e9),' Mt'])
    
%     % Relative importance of wells
%     num_init_wells = nw;
%     num_optim_wells = numel(avg_optim_masses(avg_optim_masses > 2*min(avg_optim_masses)));
%     figure;
%     bar([relative_contribution_init./num_init_wells; relative_contribution_optim./num_optim_wells]')
%     xlabel('Well', 'fontsize', 16);
%     ylabel('Contribution to total injected mass (%)', 'fontsize', 16);
%     set(gca, 'fontsize', 16); % increase fontsize on axes.
%     box off;

    %% Plot mapPlot of wells:
    cinx_inj = [schedule.control(1).W.cells];
    figure;
    mapPlot(gcf, Gt, 'traps', ta.traps, ...
        'trapcolor', [0.5 0.5 0.5], 'trapalpha', 0.7, ...
        'rivers', ta.cell_lines, 'rivercolor', [1 0 0], ...
        'maplines', 20, 'wellcells', cinx_inj, 'well_numbering', true);
    colorizeCatchmentRegions(Gt, ta);
    plotGrid(Gt, 'FaceColor','none','EdgeAlpha',0.1)
    axis equal tight off
    
    % remove well number and dot of those wells deemed to have close to
    % zero-rate
    masses_tmp = avg_optim_masses;
    masses_tmp(avg_optim_masses < 2*min(avg_optim_masses)) = 0;
    hfig = gcf;
    wellLabels = findobj(hfig.Children.Children,'Type','Text');
    lines = findobj(hfig.Children.Children,'Type','Line');
    assert(numel(lines(end).XData) == nw);
    Xdots_updated = lines(end).XData;
    Ydots_updated = lines(end).YData;
    masses_tmp_flipped = fliplr(masses_tmp); % ordering of wellLabels is backwards
    for i = 1:numel(masses_tmp) 
        if masses_tmp_flipped(i) == 0
            wellLabels(i).String='  '; 
        end
    end
    for i = 1:numel(masses_tmp) % ordering of well dots is forwards
        if masses_tmp(i) == 0
            Xdots_updated(i) = 0;
            Ydots_updated(i) = 0;
        end
    end
    lines(end).XData = Xdots_updated(Xdots_updated > 0);
    lines(end).YData = Ydots_updated(Ydots_updated > 0);
    title('Remaining wells')
    
    

end
function compareWellrates_viaWellSols(wellSols1, wellSols2, schedule, co2RefRho)


    % take only wellSols.qGr corresponding to injection period (control step 1)
    wellSols1 = wellSols1( schedule.step.control == 1 );
    wellSols2 = wellSols2( schedule.step.control == 1 );

    % compute avg qGr per well
    test1 = cell2mat(wellSols1);
    test2 = cell2mat(wellSols2);
    [~,nw] = size(test1);
    assert( numel(schedule.control(1).W) == nw )
    for i=1:nw
        qGr1(i) = mean([test1(:,i).qGr]); % m3/s
        qGr2(i) = mean([test2(:,i).qGr]);
        % or keep track of max and min, and include as spread in bar-plot
    end

    %
    avg_init_rates = qGr1;
    avg_opt_rates  = qGr2;

    val = [avg_init_rates; avg_opt_rates]';
    
    % values are currently in m3/sec using the CO2 reference density.  We
    % convert it to Mt/year
    val = val * co2RefRho * year / 1e6 / 1e3;
    
    figure; bar(val);
    set(gca, 'xlim', [0 size(val, 1) + 1]);
    xlabel('Well', 'fontsize', 16);
    ylabel('Rate (Mt/year)', 'fontsize', 16);
    set(gca, 'fontsize', 16); % increase fontsize on axes.
    box off;


end
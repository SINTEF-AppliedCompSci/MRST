function compareWellrates(initSchedule, optimisedSchedule, co2RefRho)

   init_rates = [initSchedule.control(1).W.val];
   opt_rates  = [optimisedSchedule.control(1).W.val];

   val = [init_rates; opt_rates]';
    
    % values are currently in m3/sec using the CO2 reference density.  We
    % convert it to Mt/year
    val = val * co2RefRho * year / 1e6 / 1e3;
    
    figure; bar(val);
    set(gca, 'xlim', [0 size(val, 1) + 1]);
    xlabel('Well', 'fontsize', 14);
    ylabel('Rate (Mt/year)', 'fontsize', 14);
    set(gca, 'fontsize', 14); % increase fontsize on axes.
    box off;
end



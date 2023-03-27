function locusPlot(states, colors, theta, tgrad, first_time, downwards)

    persistent p_range;
    persistent t_range;
    
    persistent CO2precise;
    CO2precise = CO2props('rho_huge','');

    % (re-)initialize ranges;
    if first_time
        p_range = []; t_range = [];
    end
    stateinfo = [states.info];
    P = [stateinfo.intPress, stateinfo.topPress];
    T = [stateinfo.intTemp,  stateinfo.topTemp];
    p_range = max_range(p_range, P);
    t_range = max_range(t_range, T);

    %@@ hack
     % p_range = [7.7 10.2]*1e6;
     % t_range = [308 314];
     p_range = [6 10.2]*1e6;
     t_range = [300 314];
    %p_range = [0.6 13]*1e6;
    %t_range = [280 350];
    
    sign = 1;
    if nargin > 5 && ~downwards
        sign = -1; % integrate upwards from reference surface
    end
    
    %plotCO2locus(P, T, CO2, p_range, t_range);
    %plotEtaIntegral(p_range, t_range, CO2, theta, tgrad, max(max([states.h])), 1);
    plotEtaIntegral(p_range, t_range, CO2precise, theta, tgrad, sign*max(max([states.h])), 5);
    % plotting locus
    hold on;
    %for i = 1:size(P, 2)
    for i = 1
        line(T(:,i), P(:,i), 'color', [.5 .5 .5], 'LineWidth', 2);
    end
    
    xlabel('Temp. (°K)')
    ylabel('Pressure (MPa)');
    
    hold off;
    
end

function plot_total_savings_function( h, C_invest, operation_rate, tax_credit_rate, alpha )
% Plot total well cost function

    if ~ishandle(h)
        figure(h)
    else
        % Avoid stealing focus if figure already exists
        set(0, 'CurrentFigure', h); clf(h);
    end
    hold on
    xlabel('CO2 injected (Mt)')
    ylabel({['total savings (million USD)'];['assuming no leakage']})
    
    
    M_inj       = 0:1e4:2e6;                                                % tonnes injected
    C_operate   = operation_rate.*M_inj;                                    % operation cost, USD
    

    % Linear total well cost:
    x = M_inj;                      % tonnes of co2 injected by a well
    y = C_invest + C_operate;       % total cost of the well
    p = tax_credit_rate .* M_inj;   % profit of the well (assuming no leakage)
    Y = p - y;                      % total savings made

    plot(x/1e6,Y/1e6)               % plot: total savings made

    % Non-linear total well cost:
    M_crit = C_invest/(tax_credit_rate - operation_rate);                   % tonne, a critical mass required to inject for the well to be a worthwhile investment
    %M_crit = C_invest/(operation_rate);
    
    alpha_max = C_invest/((tax_credit_rate-operation_rate)*M_crit);
    
    
    y = C_invest .* tanh(M_inj/(alpha*M_crit)) + C_operate;                 % total cost of the well
    Y = p - y;                                                              % total savings made
    
    plot(x/1e6,Y/1e6, '--r', 'LineWidth',3)                                 % plot: total savings made
    
    %plot((M_crit)/1e6, (C_invest*tanh(1/alpha)+operation_rate*M_crit)/1e6, '*')
    ylim = get(gca,'ylim');
    plot([M_crit M_crit]./1e6, ylim, '--')                                  % plot: critical line
    grid on

    legend('linear','non-linear','critical line')
    
end
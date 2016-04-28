function plot_total_well_cost_function( h, C_invest, operation_rate, tax_credit_rate, alpha )
% Plot total well cost function

    if ~ishandle(h)
        figure(h)
    else
        % Avoid stealing focus if figure already exists
        set(0, 'CurrentFigure', h); clf(h);
    end
    hold on
    xlabel('CO2 injected (Mt)')
    ylabel('total well cost (million USD)')
    
    
    M_inj       = 0:1e4:2e6;                                                % tonnes injected
    C_operate   = operation_rate.*M_inj;                                    % operation cost, USD

    
    % Linear total well cost:
    x = M_inj;                                                      % tonnes of co2 injected by a well
    y = C_invest + C_operate;                                       % total cost of the well

    M_crit = C_invest/(tax_credit_rate - operation_rate);           % tonne, a critical mass required to inject for the well to be a worthwhile investment
    %M_crit = C_invest/(operation_rate);
    
    plot(x/1e6,y/1e6, 'b')                                           % linear total well cost
    plot((M_crit)/1e6, (C_invest+operation_rate*M_crit)/1e6, '*b')   % critical injection mass for linear total well cost
    
    
    
    % Non-linear total well cost:
    alpha_max = C_invest/((tax_credit_rate-operation_rate)*M_crit);
    
    y = C_invest .* tanh(M_inj/(alpha*M_crit)) + C_operate;         % total cost of the well

    plot(x/1e6,y/1e6, '--r', 'LineWidth',3)                                     % non-linear total well cost
    plot((M_crit)/1e6, (C_invest*tanh(1/alpha)+operation_rate*M_crit)/1e6, '*r') % critical injection mass for non-linear total well cost
    
    
    % Well profit (assuming no leakage):
    p = tax_credit_rate .* M_inj;                                   % USD
    plot(M_inj/1e6, p/1e6, '--')                                    % ideal well profit
    
    ylim = get(gca,'ylim');
    plot([M_crit M_crit]./1e6, ylim, '--')                          % break-even injection mass (where total well cost = well profit)

    grid on
    
    legend('linear','linear critical point','non-linear','non-linear critical point','ideal well profit','critical line')

end
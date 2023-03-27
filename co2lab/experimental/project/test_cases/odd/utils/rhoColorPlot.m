function [h, pvals, tvals, pgrid, tgrid] = ...
        rhoColorPlot(CO2, p_range, t_range, p_res, t_res, critpoint_sz)

    pvals = linspace(p_range(1), p_range(2), p_res);
    tvals = linspace(t_range(1), t_range(2), t_res);
    [tgrid, pgrid] = meshgrid(tvals, pvals);
    rho = reshape(CO2.rho(pgrid(:), tgrid(:)), p_res, t_res);
    h = pcolor(tvals, pvals, rho);
    set(h, 'edgeAlpha', 0);
    view(0, 90);
    hold on;
    
    % adding liquid-vapor-boundary line    
    [pc, tc] = CO2CriticalPoint();
    lvb_t = linspace(t_range(1), tc, 100)';
    lvb_p = CO2VaporPressure(lvb_t);
    
    line(lvb_t, lvb_p, 'lineWidth', 4, 'color', [0 0 0]);
    
    % Indicating critical point
    if (nargin > 5)
        plot(tc, pc, 'kx', 'MarkerSize', critpoint_sz, 'LineWidth', 4);
    end
    %plot(tc, pc, 'o', 'MarkerSize', 10, 'LineWidth', 4);
    
end

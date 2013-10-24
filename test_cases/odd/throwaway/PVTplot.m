function h = PVTplot(Prange, Trange, CO2)
%
% Plot the CO2 PV-plane, within the specified range
%
% SYNOPSIS:
%   function h = PVTplot(Prange, Trange)
%
% PARAMETERS:
%   Prange - minimum and maximum pressure for the plot (2-component vector)
%   Trange - minimum and maximum temperature for the plot (2-component vector)
%      CO2 - (optional).  If provided, use its function 'rho(P,T)' to plot
%            density lines
%
% RETURNS:
%   h - handle to graphic
%
% EXAMPLE:
%
% SEE ALSO:
%
    h = figure;
    axis([Trange, Prange]);
    hold on;

    %% getting values for critical point
    [Pc, Tc] = CO2CriticalPoint; % CO2 critical point    
    
    %% Drawing vapor-liquid boundary
    VLB_T = linspace(Trange(1), Tc, 100);
    VLB_P = CO2VaporPressure(VLB_T);
    plot(VLB_T, VLB_P, 'k');
    
    % Drawing critical point
    plot(Tc, Pc, 'ok');

    % Drawing supercritical region
    line([Tc, Tc], [Pc, Prange(2)], 'color', 'black', 'LineStyle', '--');
    line([Tc, Trange(2)], [Pc, Pc], 'color', 'black', 'LineStyle', '--');

    if nargin > 2
        res_t = 100;
        res_p = 100;
        pvals = linspace(Prange(1), Prange(2), res_p);
        tvals = linspace(Trange(1), Trange(2), res_t);
        [t, p] = meshgrid(tvals, pvals);
        rho = reshape(CO2.rho(p(:), t(:)), res_t, []);
        
        [c, h] = contour(tvals, pvals, rho, 'k:');
        clabel(c, h,'FontSize', 16);
    end
    
    % % adding labels
    xlabel('T (degree K)', 'FontSize', 20);
    ylabel('P (MPa)', 'FontSize', 20);
    set(gca, 'FontSize', 16);
    set(gca, 'yticklabel', get(gca,'ytick')/1e6);    
    % th1 = text(290, 11e6, 'liquid');
    % th2 = text(320,  6e6, 'gas');
    % th3 = text(320,  9e6, 'supercritical');
    
end

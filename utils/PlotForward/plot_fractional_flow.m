function plot_fractional_flow(sw_cell,fw_cell)
    figure('Name','Fractional flow','NumberTitle','off')
    plot(sw_cell,fw_cell, '-bo', 'Markersize', 3)
    xlim([0 1])
    ylim([0 1])
    xlabel('Water Saturation')
    ylabel('Fractional Flow')
    grid on
end
function plot_saturation_profile(s, seff, smax, G)

    % define colors
    red = [1 0 0]; blue = [0 0 1]; white = [1 1 1]; black = [0 0 0];
    mix = @(a, b, c) (a * c + b * (100 - c)) / 100;
    trapped_color = mix(mix(red, white, 40), black, 80);;
    plume_color = mix(red, white, 40);
    residual_color = mix(mix(red, blue, 50), white, 20);
    brine_color = mix(blue, white, 30);
    
    zvals = G.cells.centroids(:,3);
    hold on;
    % current saturation, colored field
    fill([s; s * 0], [zvals; flipud(zvals)], plume_color); 
    % immobilized saturation (current - effective saturation)
    fill([s - seff; 0*s], [zvals; flipud(zvals)], residual_color);
    % brine zone, colored field
    fill([s; s * 0 + 1], [zvals; flipud(zvals)], brine_color); 
    
    % historically maximal saturation
    last_ix = find(smax>0, true, 'last');
    last_ix = min(last_ix+1, numel(value(smax)));
    plot(smax(1:last_ix), zvals(1:last_ix), 'linewidth', 2, 'color', 'r');
    set(gca, 'ydir', 'reverse')

end
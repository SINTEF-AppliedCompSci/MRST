function plotsummary(abscissa, z_top, saved_tsteps, colors, column_labels, tilt_offset)
    
    columns = size(saved_tsteps, 1);
    % stylemap for twin graphs (top/iface pressure, top/iface density)
    dashed_colors = cellfun(@(x) [x, '--'], colors, 'UniformOutput', false);
    
    sym_colors = {'o', '*', '<'}
    
    colors_twin = {dashed_colors{:}, colors{:}};

    for c = 1:columns
        
        stateinfo = [saved_tsteps(c,:).info];

        
        % Plotting height
        subplot(4,columns, (1-1) * columns + c) % first row
        
        heights = bsxfun(@plus, [saved_tsteps(c,:).h], z_top);
        max_z = max(z_top);
        min_z = min(z_top);
        polygraph([heights, z_top], {colors{:}, 'k'}, ...
                  {'km','m'}, '', abscissa/1e3, ...
                  [min_z, max_z + (max_z-min_z)* 0.1]);

        polygraph([heights(1:6:end,1)], {sym_colors{1}, 'k'}, ...
                  {'km','m'}, '', abscissa(1:6:end)/1e3, ...
                  [min_z, max_z + (max_z-min_z)* 0.1]);
        polygraph([heights(3:6:end,2)], {sym_colors{2}, 'k'}, ...
                  {'km','m'}, '', abscissa(3:6:end)/1e3, ...
                  [min_z, max_z + (max_z-min_z)* 0.1]);
        polygraph([heights(5:6:end,3)], {sym_colors{3}, 'k'}, ...
                  {'km','m'}, '', abscissa(5:6:end)/1e3, ...
                  [min_z, max_z + (max_z-min_z)* 0.1]);
        
        set(gca, 'YDir', 'reverse');

        % Plotting pressures
        subplot(4,columns, (2-1) * columns + c) % second row
        polygraph([stateinfo.topPress, stateinfo.intPress]/1e6, ...
                  colors_twin, {'km', 'MPa'}, ...
                  '', abscissa/1e3, []);
        
        % plotting densities
        subplot(4,columns, (3-1) * columns + c) % third row
        polygraph([stateinfo.topRho, stateinfo.intRho], colors_twin, {'km', 'kg/m^3'}, ...
                  '', abscissa/1e3, []);
        
        % Plotting mass flux
        subplot(4,columns, (4-1) * columns + c) % fourth row
        polygraph([stateinfo.fluxCO2], colors, {'km','kg/sm^2'}, '', abscissa(2:end)/1e3, []);

        
    end    

    % column labeling
    for c = 1:columns
        subplot(4, columns, (1-1)*columns + c); 
        xlim = get(gca, 'Xlim');
        ylim = get(gca, 'Ylim');
        xpos = xlim(1) + diff(xlim)/2;
        ypos = ylim(1) - diff(ylim)/5;
        h = text(xpos, ypos, column_labels{c}, 'FontSize', 20, 'HorizontalAlignment', 'center');
    end
    
    % row labeling
    row_labels = {'Depth', 'Pressure', 'Density', 'Mass flux'};
    for r = 1:4
        subplot(4, columns, (r-1)* columns + 1)
        xlim = get(gca, 'Xlim');
        ylim = get(gca, 'Ylim');
        xpos = xlim(1) - diff(xlim)/3;
        ypos = ylim(1) + diff(ylim)/2;
        h = text(xpos, ypos, row_labels{r}, 'FontSize', 20, 'HorizontalAlignment', ...
                 'center', 'rotation', 90);
    end
    
    % Adding second axis to height plots if required
    if exist('tilt_offset') && tilt_offset ~= 0
        for c = 1:columns
            subplot(4, columns, c);
            addSecondAxis(gca, tilt_offset);
        end
    end
end


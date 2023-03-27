function yrange = polygraph(graphs, colors, labels, plot_title, ...
                            xvals, yrange, include_in_legend, keep_axis)
    if ~exist('keep_axis')
        keep_axis = false;
    end
    if ~exist('include_in_legend')
        include_in_legend = repmat(true, size(graphs,2), 1);
    end
    
    if isempty(graphs)
        yrange = [];
        return
    end
        
    hold on; %cla;
    for g_ix = 1:size(graphs, 2) % one graph per column
        
        h = plot(xvals, graphs(:, g_ix), colors{g_ix});
        if ~include_in_legend(g_ix)
            legendOff(h);
        end
    end
    xlabel(labels{1}, 'fontsize', 15);
    ylabel(labels{2}, 'fontsize', 15);

    if isempty(yrange) || max(max(graphs)) > yrange(2) || min(min(graphs)) < yrange(1)
        % Y-range not yet set/imposed.  Computing a reasonable range
        PADDING = 0.1;
        
        yrange = [min(min(graphs)), max(max(graphs))];
        span = yrange(2) - yrange(1);
        yrange(1) = yrange(1) - PADDING * span;
        yrange(2) = yrange(2) + PADDING * span;
    
    end        
   
    title(plot_title);
    
    % small 'hack' to avoid collapse if y is constant (range of zero...)
    if yrange(1) < yrange(2)
        margin = 0;
    else
        if abs(yrange(1)) > 100 * eps;
            margin = yrange(1) * 0.05;
        else
            margin = 0.01;
        end
    end
    if ~keep_axis
        axis([min(xvals), max(xvals), yrange(1)-margin, yrange(2)+margin]);
    end
end

function legendOff(h)

    set(get(get(h, 'Annotation'), 'LegendInformation'), ...
        'IconDisplayStyle', 'off');
end

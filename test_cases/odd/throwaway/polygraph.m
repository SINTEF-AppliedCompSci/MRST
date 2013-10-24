function yrange = polygraph(graphs, colors, labels, plot_title, xvals, yrange)

    hold on; cla;
    for g_ix = 1:size(graphs, 2) % one graph per column
        
        plot(xvals, graphs(:, g_ix), colors{g_ix});
        
    end
    xlabel(labels{1});
    ylabel(labels{2});

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
        margin = yrange(1) * 0.05;
    end
    axis([min(xvals), max(xvals), yrange(1)-margin, yrange(2)+margin]);
    
end

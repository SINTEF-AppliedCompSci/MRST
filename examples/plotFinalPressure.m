function plotFinalPressure(G, states, name)

    figure
    plotToolbar(G, states)
    axis tight equal
    colorbar
    title(name)

    figure, hold on
    plotCellData(G, states{end}.pressure, 'edgealpha', 0);
    contour(reshape(G.cells.centroids(:, 1), G.cartDims), ...
            reshape(G.cells.centroids(:, 2), G.cartDims), ...
            reshape(states{end}.pressure, G.cartDims), ...
            'linewidth', 1, 'color', 'k');
    axis tight equal
    colorbar
    title([name, ' at endtime'])

end

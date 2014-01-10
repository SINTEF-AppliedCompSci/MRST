function plotWellsPrint(G, W, D)
% Plots wells as simple colored circles. Helper for diagnostics examples.
    gc = G.cells.centroids;
    hold on
    for i = 1:numel(W)
        c = W(i).cells(1);
        if ~ismember(i, D.inj)
            color = 'red';
        else
            color = 'blue';
        end
        plot3(gc(c, 1), gc(c, 2), -5, 'O', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', color, 'MarkerSize', 12)
    end
end

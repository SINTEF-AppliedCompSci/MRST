function plotGmsh(G)

    rng(1);

    if G.griddim == 2

        % Cells
        figure, hold on
        ut = unique(G.cells.tag);
        colors = [lines(7); rand(numel(ut)-7, 3)];

        for k = 1:numel(ut)
            idx = G.cells.tag == ut(k);
            plotGrid(G, idx, 'facecolor', colors(k, :));
        end

        % Faces
        figure, hold on
        plotGrid(G)
        ut = unique(G.faces.tag);
        colors = [lines(7); rand(numel(ut)-7, 3)];

        for k = 1:numel(ut)
            idx = G.faces.tag == ut(k);
            plotFaces(G, idx, 'linewidth', min(2, k), 'edgecolor', colors(k, :));
        end

    elseif G.griddim == 3

        % Cells
        figure; hold on
        ut = unique(G.cells.tag);
        colors = [lines(7); rand(numel(ut)-7, 3)];
        plotGrid(G, 'facecolor', 'none', 'facealpha', 0.0)
        view(3)

        for k = 1:numel(ut)
            if ut(k) > 0
                idx = G.cells.tag == ut(k);
                plotGrid(G, idx, 'facecolor', colors(k, :));
            end
        end

        % Faces
        fig = figure; hold on;
        ut = unique(G.faces.tag);
        colors = [lines(7); rand(numel(ut)-7, 3)];
        plotGrid(G, 'facecolor', 'none', 'facealpha', 0.0)

        for k = 1:numel(ut)
            if ut(k) > 0
                figure
                idx = G.faces.tag == ut(k);
                plotFaces(G, idx, 'linewidth', 2, 'facecolor', colors(k, :));
                view(3)
                title(sprintf('tag %d: %d entities', ut(k), sum(idx)));

                figure(fig)
                plotFaces(G, idx, 'linewidth', 2, 'facecolor', colors(k, :));
            end
        end

    end

end
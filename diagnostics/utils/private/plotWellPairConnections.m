function varargout = plotWellPairConnections(G, WP, D, W, pv)
    hold on

    maxAlloc = 0;
    for i = 1:numel(WP.inj)
        maxAlloc = max([maxAlloc sum(WP.inj(i).alloc, 1)]);
    end

    handles = [];
    for i=1:numel(D.inj)
        ipos = mean(G.cells.centroids(W(D.inj(i)).cells,:), 1);
        ialloc = sum(WP.inj(i).alloc, 1);


        for p=1:numel(D.prod)

            ppos = mean(G.cells.centroids(W(D.prod(p)).cells,:), 1);

            alloc = ialloc(p);


            if alloc/sum(ialloc) > 0.01

                localRegion = find(D.ipart == i & D.ppart == p);
                center = sum(G.cells.centroids(localRegion, :).*repmat(pv(localRegion), 1, G.griddim), 1)/sum(pv(localRegion));
                pts = [ipos; center; ppos];

                htext = text(center(1), center(2), center(3), ...
                             sprintf('%2.1f%%', 100*(alloc/sum(ialloc))),...
                             'VerticalAlignment', 'Bottom', 'FontWeight', 'bold', 'FontSize', 12, 'Color', colorizeWell('inj', i, D));

                thick = 20*alloc/(maxAlloc);

                hline = plot3(pts(:,1), pts(:,2), pts(:,3),...
                       '-','LineWidth', thick, 'Color', colorizeWell('prod', p, D));
                handles = [handles; hline; htext]; %ok
            end
        end
    end

    if nargout > 0
        varargout{1} = handles;
    end
    hold off
end

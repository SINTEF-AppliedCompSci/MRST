function wcells = findQuarter5spotLocations(G,p1,p2)
d1 = sqrt((G.cells.centroids(:,2)-p1(2)).^2+(G.cells.centroids(:,1)-p1(1)).^2);
d2 = sqrt((G.cells.centroids(:,2)-p2(2)).^2+(G.cells.centroids(:,1)-p2(1)).^2);
[~,wcells(1)] = min(d1);
[~,wcells(2)] = min(d2);
end
function [pts, removed] = surfaceSufCondFromGrid3D(pts, grids_2, gamma)

removed = false(size(pts, 1), 1);
for i = 1:numel(grids_2)
   g = grids_2{i};
   g = mrstGridWithFullMappings(g);
   cells = g.nodes.cells(g.nodes.cellPos(1:end-1));
   kappaSqr = sum((g.nodes.coords - g.cells.sites(cells, :)).^2,2);
   R = sqrt(kappaSqr + gamma^2);
   [~, loc_removed] = surfaceSufCond3D(pts, g.nodes.coords, R);
   removed = removed | loc_removed;
end
pts = pts(~removed, :);

end

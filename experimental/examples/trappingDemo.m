moduleCheck co2lab

coarsening_level = 1; % grid downsampling factor (1 = no downsampling)
Gt = getFormationTopGrid('Krossfjordfm', coarsening_level);

% Identify boundary faces
bfaces = prod(Gt.faces.neighbors, 2) == 0;

% identify faces east of 5.7e5
efaces = Gt.faces.centroids(:,1) > 5.7e5;

% identify boundary faces east of 5.7e5 (these will be closed
cfaces = bfaces & efaces;

ta = trapAnalysis(Gt, false, 'closed_boundary_edges', find(cfaces));


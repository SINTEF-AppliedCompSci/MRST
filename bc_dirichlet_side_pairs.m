function bc_std = bc_dirichlet_side_pairs(G, dim, dirval)

% Set dirichlet values on sides dim. Homogeneous Neumann otherwise.
dirichlet_val_1 = dirval;
dirichlet_val_2 = -dirval;

bf = boundaryFaces(G);
bc_std = [];
isbdry = (1:G.faces.num)';
tol = 1e-5;

xmin = min(G.nodes.coords(:,dim));
xmax = max(G.nodes.coords(:,dim));
imin = G.faces.centroids(:,dim) < xmin+tol;
imax = G.faces.centroids(:,dim) > xmax-tol;
dirichlet_ix_1 = isbdry(imin);
dirichlet_ix_2 = isbdry(imax);

bc_std = addBC(bc_std, dirichlet_ix_1, 'pressure', dirichlet_val_1);
bc_std = addBC(bc_std, dirichlet_ix_2, 'pressure', dirichlet_val_2);

% Set the rest as homogeneous Neumann
bc_std = addBC(bc_std, setdiff(bf, [dirichlet_ix_1, dirichlet_ix_2]), ...
               'flux', 0.0);

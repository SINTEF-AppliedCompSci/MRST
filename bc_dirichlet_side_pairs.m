function bc_std = bc_dirichlet_side_pairs(G, dim)

% Set dirichlet values on sides dim. Homogeneous Neumann otherwise.
% NB: we want to avoid zero pressure in ntpfa:
% 
dirichlet_val_min = 1;
dirichlet_val_max = 10;

bf = boundaryFaces(G);
bc_std = [];
isbdry = (1:G.faces.num)';
tol = 1e-5;

xmin = min(G.nodes.coords(:,dim));
xmax = max(G.nodes.coords(:,dim));
imin = G.faces.centroids(:,dim) < xmin+tol;
imax = G.faces.centroids(:,dim) > xmax-tol;
dirichlet_ix_min = isbdry(imin);
dirichlet_ix_max = isbdry(imax);

bc_std = addBC(bc_std, dirichlet_ix_min, 'pressure', dirichlet_val_min);
bc_std = addBC(bc_std, dirichlet_ix_max, 'pressure', dirichlet_val_max);

% Set the rest as homogeneous Neumann
bc_std = addBC(bc_std, setdiff(bf, [dirichlet_ix_min, dirichlet_ix_max]), ...
               'flux', 0.0, 'sat', [1, 0]);

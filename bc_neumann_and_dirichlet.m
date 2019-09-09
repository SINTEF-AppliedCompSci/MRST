function bc_std = bc_neumann_and_dirichlet(G)

% Case with Neumann on xmin and Dirichlet on xmax. Hom Neumann
% otherwise.
neumann_val = 1; % Hard coded in convertBC
dirichlet_val = 0; % Hard coded in convertBC

% Don't use fluxside / pside but set manually to get correct
% scaling. 
bf = boundaryFaces(G);
bc_std = [];
isbdry = (1:G.faces.num)';
tol = 1e-5;
xmin = min(G.nodes.coords(:,1));
xmax = max(G.nodes.coords(:,1));
imin = G.faces.centroids(:,1) < xmin+tol;
imax = G.faces.centroids(:,1) > xmax-tol;
neumann_ix = isbdry(imin);
area = G.faces.areas(neumann_ix);
bc_std = addBC(bc_std, neumann_ix, 'flux', neumann_val.*area, 'sat', []);
dirichlet_ix = isbdry(imax);
bc_std = addBC(bc_std, dirichlet_ix, 'pressure', dirichlet_val, ...
               'sat', []);
% Set the rest as hom Neumann
bc_std = addBC(bc_std, setdiff(bf, union(neumann_ix, dirichlet_ix)), ...
               'flux', 0.0, 'sat', []);


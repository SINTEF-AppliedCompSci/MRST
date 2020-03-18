function bc_std = bc_neumann_and_dirichlet(G)

% Case with Neumann on xmin and Dirichlet on xmax. Hom Neumann
% otherwise.
neumann_val = 1; % Hard coded in convertBC
dirichlet_val = 1; % Hard coded in convertBC
sat = 1;

% Don't use fluxside / pside but set manually to get correct
% scaling.
bf = boundaryFaces(G);
bc_std = [];

isbdry = (1:G.faces.num)';
tol = 1e-5;

xmax = max(G.nodes.coords(:, 1));
imax = G.faces.centroids(:, 1) > xmax - tol;
neumann_ix = isbdry(imax);
area = G.faces.areas(neumann_ix);
bc_std = addBC(bc_std, neumann_ix, ...
    'flux', neumann_val, ...
    'sat', sat*ones(numel(neumann_val), 1));

xmin = min(G.nodes.coords(:, 1));
imin = G.faces.centroids(:, 1) < xmin + tol;
dirichlet_ix = isbdry(imin);
bc_std = addBC(bc_std, dirichlet_ix, ...
    'pressure', dirichlet_val, ...
    'sat', sat*ones(numel(dirichlet_val), 1));

% Set the rest as hom Neumann
ix = setdiff(bf, union(neumann_ix, dirichlet_ix));
bc_std = addBC(bc_std, ix, ...
    'flux', 0.0, ...
    'sat', sat*ones(numel(ix), 1));

% Set saturation on in inflow faces
%bc_std.sat = sat*ones(numel(dirichlet_ix), 1);
%bc_std.sat = sat*ones(numel(neumann_ix), 1);
%bc_std.sat = sat*ones(G.cells.num, 1);
%bc_std.sat = sat*ones(numel), 1);


% % Test
% bc = [];
% bc = pside(bc, G, 'XMax', dirichlet_val);
% bc = fluxside(bc, G, 'XMin', neumann_val);
end

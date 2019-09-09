function bc_std = bc_hom_neumann(G)

% Case with hom Neumann everywhere
bf = boundaryFaces(G);
isbdry = (1:G.faces.num)';
neumann_ix = isbdry(bf);
bc_std = [];
bc_std = addBC(bc_std, neumann_ix, 'flux', 0.0, 'sat', []);


function stable = isStable(sol, g)

%{
#COPYRIGHT#
%}


tol = 1e-6;

stable = false;

magn = @(v)(sqrt(sum(v.^2,2)));

n = g.cells.normals(:,3)./magn(g.cells.normals);

h_z = n(:,1).*sol.h;

ix = (h_z < (g.cells.H-g.cells.H*tol)) & (h_z > g.cells.H*tol);

level = h_z + g.cells.z + abs(min(g.cells.z));

u = unique(level(ix));


if (max(u) < min(u) + max(u)*tol)
   stable = true;

end

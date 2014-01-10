cartDims = [500, 500];
physDims = [1000, 1000];
g = computeGeometry(cartGrid(cartDims, physDims));

rock.perm = ones(g.cells.num, 1);

BI = mex_ip_simple(g, rock);

%{
W = addWell([], g, rock, 1);
W = addWell(W , g, rock, g.cells.num);

[BI, connPos, conns] = mex_ip_simple(g, rock, W);

nconn = diff(connPos);

BI([nconn(1)^2, end]) = [ W.WI ];

[S, r, F, L] = mex_schur_comp_symm(BI, connPos, conns);

[i, j] = blockDiagIndex(nconn, nconn);
SS = sparse(double(conns(i)), double(conns(j)), S);
R  = accumarray(conns, r);

R(end-1 : end) = 1000 * [1, -1];

lam = SS \ R;

[flux, press] = mex_compute_press_flux(BI, lam, connPos, conns, F, L);

plotCellData(g, press);
%}

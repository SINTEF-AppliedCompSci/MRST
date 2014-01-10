%run ../../startup
G = computeGeometry(cartGrid([200e3, 1], [1, 1]));
rock.perm = ones(G.cells.num, 1);

BI = mex_ip_simple(G, rock);

connPos = G.cells.facePos;
conns   = G.cells.faces(:,1);

[S, r, F, L, q] = mex_schur_comp_symm(BI, connPos, conns, []);

nconn  = diff(connPos);
[i, j] = blockDiagIndex(nconn, nconn);

SS = sparse(double(conns(i)), double(conns(j)), S);
R  = accumarray(conns, r);

SS(1) = SS(1) * 2;
R([1, G.cells.num+1]) = [1, -1];

x = SS \ R;
%%
[v, p] = mex_compute_press_flux(BI, x, connPos, conns, F, L, q);

%%
plotCellData(G, p);

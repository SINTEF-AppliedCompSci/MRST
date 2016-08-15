
G = cartGrid([5,5]);

G = computeVEM2DGeometry(G);
G = sortEdges(G);

state = initState(G, [], 0);
rock.perm = ones(G.cells.num,1);
fluid = initSingleFluid('mu', 1, 'rho', 1);


S = computeVirtualIP(G, rock, fluid, 1);
A = S.A;

ii = 1:numel(S.dofVec); jj = S.dofVec;
P = sparse(ii,jj,1, numel(S.dofVec), G.nodes.num);

A = P'*A*P;



U = A\ones(G.nodes.num,1);
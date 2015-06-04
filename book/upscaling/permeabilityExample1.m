%% Upscale a problem with a diagonal trend
% In this example we will upscale a problem with diagonal trend and compare
% the results obtained with a flow-based method with pressure-drop and
% periodic boundary conditions

mrstModule add incomp upscaling;

%% Set up model
[Lx,Ly] = deal(10);
[nx,ny] = deal(20);
G = computeGeometry(cartGrid([nx ny], [Lx Ly]));
permX = 30*milli*darcy;
[i, j] = gridLogicalIndices(G);
poro = 0.3;
rock.perm = repmat(30*milli*darcy, [G.cells.num, 1]);
rock.poro = repmat(.3, [G.cells.num, 1]);
rock.perm( sin(2*sum(G.cells.centroids(:,:),2)*pi/Lx)>0 ) = 10*milli*darcy;
hT    = computeTrans(G, rock);
fluid = initSingleFluid('mu' ,1, 'rho', 1);

clf, plotCellData(G,rock.perm/(milli*darcy),'EdgeColor','none'); colorbar

%% Structures with boundary conditions
d = G.griddim;
[bcl,bcr, Dp]=deal(cell(d,1));
bcsides = {'XMin', 'XMax'; 'YMin', 'YMax'; 'ZMin', 'ZMax'};
for j = 1:d;
   bcl{j} = pside([], G, bcsides{j, 1}, 0);
   bcr{j} = pside([], G, bcsides{j, 2}, 0);
   Dp{j}  = 0;
end
Dp{1} = 4*barsa;
L  = max(G.faces.centroids)-min(G.faces.centroids);


%% Pressure-drop boundary conditions
[v,dp] = deal(zeros(d, 1));
for i=1:d
   bc = addBC([], bcl{i}.face, 'pressure', Dp{1});
   bc = addBC(bc, bcr{i}.face, 'pressure', Dp{2});

   xr = initResSol(G, 100*barsa, 1);
   xr = incompTPFA(xr, G, hT, fluid, 'bc', bc);

   v(i)  = sum(xr.flux(bcr{i}.face)) / ...
           sum(G.faces.areas(bcr{i}.face));
   dp(i) = Dp{1}/L(i);
end
K = convertTo(v./dp, milli*darcy)                                          %#ok<NASGU,NOPTS>


%% Periodic boundary conditions
[Gp, bcp] = makePeriodicGridMulti3d(G, bcl, bcr, Dp);
ofaces = cell(d,1);
for j=1:d, ofaces{j} = bcp.face(bcp.tags==j); end
v  = nan(d);
dp = Dp{1}*eye(d);
nbcp = bcp;
for i=1:d
   for j=1:d, nbcp.value(bcp.tags==j) = dp(j,i); end

   xr = initResSol(Gp, 100*barsa, 1);
   xr = incompTPFA(xr, Gp, hT, fluid, 'bcp', nbcp);

   for j=1:d
      v(j,i) = sum(xr.flux(ofaces{j})) / ... *bcp.sign(bcp.tags==j)) / ...
               sum(Gp.faces.areas(ofaces{j}));
   end
end
dp = bsxfun(@rdivide, dp, L);
K = convertTo(v/dp, milli*darcy)                                          %#ok<NOPTS>
[V,D] = eig(K)                                                             %#ok<NOPTS>
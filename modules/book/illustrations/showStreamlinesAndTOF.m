mrstModule add diagnostics incomp streamlines;

%% Set up and solve flow problem
[nx,ny] = deal(64);
G = cartGrid([nx,ny],[500,250]);
G = computeGeometry(G);
p = gaussianField(G.cartDims, [0.2 0.4], [11 3], 2.5);
K = p.^3.*(1.5e-5)^2./(0.81*72*(1-p).^2);
rock.poro = p(:);
rock.perm = K(:);
hT = computeTrans(G, rock);

%%
gravity reset off
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);
pv  = sum(poreVolume(G,rock));

n = 12;
src = addSource([],  nx*n+n+1, pv/2);
src = addSource(src, nx*n+n+1+nx-2*n, pv/2);
src = addSource(src, G.cells.num-(n-.5)*nx, -pv);
state = initResSol(G, 0.0, 1.0);
state = incompTPFA(state, G, hT, fluid, 'src', src);


%% Compute time-of-flight
clf
tof = computeTimeOfFlight(state, G, rock, 'src', src);
x = reshape(G.cells.centroids(:,1),G.cartDims);
y = reshape(G.cells.centroids(:,2),G.cartDims);
z = reshape(tof,G.cartDims);
contour(x,y,z,exp([-4.5:.5:-1.5 -1.25:.25:0.25 .5:.25:2.5]),'Color',[.6 .6 .6]);
axis tight; box on; set(gca,'XTick',[],'YTick',[]);

%% Trace streamlines
seed = (nx*ny/2 + (1:2:nx)).';
Sf = pollock(G, state, seed, 'substeps', 1);
Sb = pollock(G, state, seed, 'substeps', 1, 'reverse', true);
hf=streamline(Sf);
hb=streamline(Sb);
set([hf; hb],'Color','k','LineWidth',1.5);

%% Add local block
clf
flag = false(G.cartDims);
flag(28:36,29:41)=true;
plotCellData(G,.3*rock.perm(:,1), 'EdgeColor','none');
plotCellData(G,rock.perm(:,1), flag(:), 'EdgeColor','none');
colormap(.7*flipud(gray)+.3*ones(size(gray)));
hf=streamline(Sf);
hb=streamline(Sb);
set([hf; hb],'Color','k','LineWidth',1);
n = gridCellNodes(G,find(flag(:)));
m = min(G.nodes.coords(n,:));
M = max(G.nodes.coords(n,:));
hold on,
plot([m(1) M([1 1]) m([1 1])], [m([2 2]) M([2 2]) m(2)],'k--','LineWidth',1);
axis tight off
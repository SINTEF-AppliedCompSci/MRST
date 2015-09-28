mrstModule add diagnostics incomp streamlines;

%% Set up and solve flow problem
[nx,ny] = deal(64);
G = cartGrid([nx,ny,1],[500,250,10]);
G = computeGeometry(G);
p = gaussianField(G.cartDims(1:2), [0.2 0.4], [11 3], 2.5);
K = p.^3.*(1.5e-5)^2./(0.81*72*(1-p).^2);
rock.poro = p(:);
rock.perm = K(:);
hT = computeTrans(G, rock);

%% Set up and solve flow problem, compute diagnostics
gravity reset off
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);
pv  = sum(poreVolume(G,rock));

n = 12;
W = addWell([],  G, rock, nx*n+n+1, ...
    'Type', 'rate', 'Comp_i', 1, 'name', 'I1', 'Val', pv/2);
W = addWell(W, G, rock, nx*n+n+1+nx-2*n, ...
    'Type','rate',  'Comp_i', 1, 'name', 'I2', 'Val', pv/2);
W = addWell(W, G, rock, round(G.cells.num-(n/2-.5)*nx-.3*nx), ...
    'Type','rate',  'Comp_i', 0, 'name', 'P1', 'Val', -pv/4);
W = addWell(W, G, rock, G.cells.num-(n-.5)*nx, ...
    'Type','rate',  'Comp_i', 0, 'name', 'P2', 'Val', -pv/2);
W = addWell(W, G, rock, round(G.cells.num-(n/2-.5)*nx+.3*nx), ...
    'Type','rate',  'Comp_i', 0, 'name', 'P3', 'Val', -pv/4);
state = initState(G, W, 0.0, 1.0);
state = incompTPFA(state, G, hT, fluid, 'wells', W);
D = computeTOFandTracer(state, G, rock, 'wells', W);

%% Trace streamlines
seed = (nx*ny/2 + (1:nx)).';
Sf = pollock(G, state, seed, 'substeps', 1);
Sb = pollock(G, state, seed, 'substeps', 1, 'reverse', true);

%% Forward time-of-flight
clf
hf=streamline(Sf);
hb=streamline(Sb);
plotGrid(G,vertcat(W.cells),'FaceColor','none','EdgeColor','r','LineWidth',1.5);
plotWell(G,W,'FontSize',20); 
set([hf; hb],'Color','w','LineWidth',1.5);
hd=plotCellData(G,D.tof(:,1),'EdgeColor','none','FaceAlpha',.6);
colormap(parula(32));
axis equal off;
print -dpng diagnost-ftof.png;

%% Backward time-of-flight
delete(hd);
hd=plotCellData(G,D.tof(:,2),'EdgeColor','none','FaceAlpha',.6);
print -dpng diagnost-btof.png;

%% Total travel time/residence time
delete(hd);
hd=plotCellData(G,sum(D.tof,2),'EdgeColor','none','FaceAlpha',.6);
print -dpng diagnost-ttof.png;

%% Tracer from I1
delete(hd);
t = D.itracer(:,1);
hd = plotCellData(G, t, t>1e-6, 'EdgeColor','none','FaceAlpha',.6);
print -dpng diagnost-C-I1.png;

%% Tracer from P2
delete(hd);
t = D.ptracer(:,2);
hd = plotCellData(G, t, t>1e-6, 'EdgeColor','none','FaceAlpha',.6);
print -dpng diagnost-C-P2.png;

%% Flooded regions
delete(hd);
hd = plotCellData(G,D.ipart,'EdgeColor','none','FaceAlpha',.6);
print -dpng diagnost-ipart.png;

%% Drainage regions
delete(hd);
hd = plotCellData(G,D.ppart,'EdgeColor','none','FaceAlpha',.6);
print -dpng diagnost-ppart.png;

%%
delete(hd);
hd = plotCellData(G,D.ipart, D.ppart~=2, 'EdgeColor','none','FaceAlpha',.6);
print -dpng diagnost-wreg.png;
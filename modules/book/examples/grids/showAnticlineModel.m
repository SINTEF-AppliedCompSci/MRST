%% Simple Model of an Anticline
% Conceptual illustration of a reservoir model with aggressive coarsening
% in the aquifer zone, modest coarsening above the initial water contact,
% and original resolution for cells with low residence time

%%
mrstModule add coarsegrid diagnostics incomp spe10

%% Cartesian grid and rock parameters
[nx,ny,nz] = deal(60,60,15);
G = cartGrid([nx ny nz], [nx ny nz].*[20 10 2]*ft);

rock = getSPE10rock(1:nx,1:ny,nz:-1:1);
rock.poro(rock.poro==0) = 1e-5;

%% Make anticline structure
x = G.nodes.coords(:,1);
y = G.nodes.coords(:,2);

normalize = @(x) (x-min(x))/(max(x)-min(x));
x = 2*normalize(x)-1;
y = 2*normalize(y)-1;
r = min(sqrt(x.^2 + y.^2),1);
dz = r.^2 ./ (r.^2 + (1-r).^4);

G.nodes.coords(:,3) = G.nodes.coords(:,3) + dz*(nz+1)*2*ft;
G = computeGeometry(G);
%%
clf, set(gcf,'Position', [230 190 1200 630]);
subplot(1,3,1)
plotCellData(G,log10(rock.perm(:,1)),'EdgeColor','none');
view(-54,30); set(gca,'dataasp',[1 1 .4]), axis tight off

%% Set initial saturation
s = zeros(G.cells.num,1);
s(G.cells.centroids(:,3)>(nz+1)*2*ft) = 1;
%[nodes,pos] = gridCellNodes(G, (1:G.cells.num).');
%c = rldecode(1:G.cells.num, diff(pos), 2).';
%A = sparse(c,nodes,1)*[s, ones([G.nodes.num,1])];
%s = bsxfun(@rdivide, A(:,1:(end-1)), A(:,end));
subplot(1,3,2), cla
plotCellData(G,1-s,'EdgeAlpha',.1); caxis([0 1]);
view(-54,30); set(gca,'dataasp',[1 1 .4]), axis tight off

%% Set wells
% producers
args = {'Type', 'bhp', 'Val', 100*barsa, 'Comp_i', [0 1], 'sign', -1};
W = verticalWell([], G, rock, 28, 30, [], args{:}, 'name', 'P1');
W = verticalWell(W,  G, rock, 34, 32, [], args{:}, 'name', 'P2');

% injectors
args = {'Type', 'bhp', 'Val', 150*barsa, 'Comp_i', [1 0], 'sign', 1};
W = verticalWell(W,  G, rock, 18, 18, [], args{:}, 'name', 'I1');
W = verticalWell(W,  G, rock, 18, 44, [], args{:}, 'name', 'I2');
W = verticalWell(W,  G, rock, 44, 44, [], args{:}, 'name', 'I4');
subplot(1,3,1)
plotWell(G,W,'Color','k','FontSize',12); axis tight off

%% Incompressible pressure solution and time-of-flight
hT = computeTrans(G, rock);
fluid = initSimpleFluid('mu' , [1 1]*centi*poise     , ...
                        'rho', [1014, 859]*kilogram/meter^3, ...
                        'n'  , [ 2, 2]);
state = initState(G, W, 100*barsa, [s 1-s]);
state = incompTPFA(state, G, hT, fluid, 'wells', W);

Tf = computeTimeOfFlight(state, G, rock, 'wells', W, 'reverse', false);
Tb = computeTimeOfFlight(state, G, rock, 'wells', W, 'reverse', true);
T  = Tf + Tb;

%% Extract the high-flow zones
% Selection criterion: travel time less than two years (times magic factor)
tfac = median(T(vertcat(W.cells)));
I = T < 2*tfac;
subplot(1,3,3), cla
plotCellData(G,I+0,'EdgeAlpha',.1); caxis([-3 2]);
plotWell(G,W,'Color','w','FontSize',12);
view(-54,30); set(gca,'dataasp',[1 1 .4]), axis tight off

%% Coarsen model differently in aquifer, oil, and high-flow zones
% Aquifer zone
pc = partitionCartGrid([nx,ny,nz],[10 10 5]);
p = compressPartition(pc);

% Zone above intial water contact
pf = partitionCartGrid([nx,ny,nz],[nx,ny,nz]/2);
p(s<1e-5) = pf(s<1e-5);

% Zone with low residence time
po = 1:G.cells.num;
p(I) = po(I);

% Near-well region
q = zeros(G.cells.num,1);
q(vertcat(W.cells))=1;
qc = accumarray(pc,q);
ind = qc(pc)>0;
p(ind) = po(ind);
subplot(1,3,1)
h = outlineCoarseGrid(G, p,'FaceColor','none','EdgeColor','k');

%%
for i=1:3
    subplot(1,3,i),
    set(gca,'Clipping', 'off', ...
        'Position',get(gca,'Position')+[-.025 -.025 .05 .05]); zoom(1.3*1.1);
end
colormap(.8*jet(128)+.2*ones(128,3));

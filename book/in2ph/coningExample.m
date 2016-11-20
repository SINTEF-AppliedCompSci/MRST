%% Simulate model with water coning
mrstModule add incomp coarsegrid


%% Make grid and assign petrophysical properties
G = cartGrid([60,40,10],[1500 1000 200]);
G.nodes.coords(:,3) = G.nodes.coords(:,3)+2050;
G = computeGeometry(G);

[x0,x1,z0,z1] = deal(675,1250,2050,2250);
flt = @(c) (c(:,1)-x0)*(z1-z0)/(x1-x0) + z0 - c(:,3);

rock.poro = .2*ones(G.cells.num,1);
rock.perm = ones(G.cells.num,1)*100*milli*darcy;
rock.perm(flt(G.cells.centroids)>0) = 10*milli*darcy;
rock.poro(flt(G.cells.centroids)>0) = 0.1;

hT = computeTrans(G,rock);

pargs = {'EdgeAlpha',.1,'EdgeColor','k'};
clf, hs = plotCellData(G,rock.perm,pargs{:});
view(3), axis tight

%% Setup wells
x = G.cells.centroids(:,[1 3]);
W = addWell([], G, rock, find(sum(bsxfun(@minus,x,[67.5 2060]).^2,2)<320), ...
            'InnerProduct', 'ip_tpf', ...
            'Type', 'bhp', 'Val', 100*barsa, ...
            'Comp_i', [0 1], 'Name', 'P', 'Dir','y');

x = G.cells.centroids(:,1:2);
W = addWell(W,  G, rock, find(sum(bsxfun(@minus,x,[1437.5 487.5]).^2,2)<320),  ...
            'InnerProduct', 'ip_tpf',...
            'Type', 'bhp', 'Val', 700*barsa, ...
            'Comp_i', [1 0], 'Name', 'I', 'Dir','z');

plotWell(G,W,'height',50,'radius',.01);
%plotGrid(G,W(1).cells,'FaceColor','g');
%plotGrid(G,W(2).cells,'FaceColor','b');

CG = generateCoarseGrid(G,(flt(G.cells.centroids)>0)+1);
plotFaces(CG,1:CG.faces.num,'FaceColor','none','LineWidth',1);
plotFaces(CG,11,'FaceColor','y','FaceAlpha',.3);
set(hs,'FaceAlpha',.35);
zoom(1.4); set(gca,'dataasp',[2 2 1]); view(25,30);

%% Fluid model
fluid = initSimpleFluid('mu' , [   1,  10] .* centi*poise     , ...
                        'rho', [1000, 100] .* kilogram/meter^3, ...
                        'n'  , [   2,   2]);

%% Simulation loop
N  = 400;
T  = 20000*day();
dT = T/N*ones(N,1);
dT = [dT(1)*sort(2.^-[1:4 4])'; dT(2:end)];

gravity reset on
rSol = initState(G, W, 0, [0, 1]);
rSol = incompTPFA(rSol, G, hT, fluid, 'wells', W);

t = 0; 
colormap(flipud(winter))
wellSols = cell(numel(dT),1);
%
for i=1:numel(dT),
   rSol = implicitTransport(rSol, G, dT(i), rock, fluid, 'wells', W);

   % Check for inconsistent saturations
   assert(max(rSol.s(:,1)) < 1+eps && min(rSol.s(:,1)) > -eps);

   % Update solution of pressure equation.
   rSol  = incompTPFA(rSol , G, hT, fluid, 'wells', W);

   % Measure water saturation in production cells in saturation
   wellSols{i} = getWellSol(W, rSol, fluid);

   % Increase time
   t = t + dT(i);

   % Plot saturation
   delete(hs);
   hs = plotCellData(G, rSol.s(:,1), (rSol.s(:,1)>.01), pargs{:});
   title([num2str(convertTo(t,day)),  ' days']),
   caxis([0 1]); drawnow
end

figure,
p = cellfun(@(x) x(1).qO', wellSols,'UniformOutput',false);
p = vertcat(p{:});
[Ym,Tm]=meshgrid(G.cells.centroids(W(1).cells,2),cumsum(dT));
plot3(Ym,Tm,p,'-ok','MarkerSize',5,'MarkerFaceColor',[.5 .5 .5]);
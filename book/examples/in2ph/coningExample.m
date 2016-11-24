%% Simulate model with water coning
mrstModule add incomp coarsegrid ad-core


%% Make grid and assign petrophysical properties
G = cartGrid([60,40,10],[1500 1000 200]);
G.nodes.coords(:,3) = G.nodes.coords(:,3)+2050;
G = computeGeometry(G);

[x0,x1,z0,z1] = deal(675,1250,2050,2250);
flt = @(c) (c(:,1)-x0)*(z1-z0)/(x1-x0) + z0 - c(:,3);

rock.poro = .2*ones(G.cells.num,1);
rock.perm = ones(G.cells.num,1)*500*milli*darcy;
rock.perm(flt(G.cells.centroids)>0) = 50*milli*darcy;
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
N  = 450;
T  = 4500*day();
dT = T/N*ones(N,1);
dT = [dT(1)*sort(2.^-[1:4 4])'; dT(2:end)];

gravity reset on
rSol = initState(G, W, 0, [0, 1]);
rSol = incompTPFA(rSol, G, hT, fluid, 'wells', W);

t = 0; 
colormap(flipud(winter))
wellSols = cell(numel(dT),1);
set(gca,'XTick',[],'Ytick',[],'ZTick',[]);
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
   % print('-dpng','-r0',sprintf('img/con-%03d.png',i));
end

%% Oil rate, with peak production indicated by red line
figure,
[Ym,Tm] = meshgrid(G.cells.centroids(W(1).cells,2),cumsum(dT)/year);
p       = cellfun(@(x) abs(x(1).qO)', wellSols,'UniformOutput',false);
p       = vertcat(p{:})*day;
[~,j]   = max(p);
m       = sub2ind(size(p),j,1:numel(W(1).cells));
surf(Ym,Tm,p); shading interp; colormap(parula(20));
hold on; plot3(Ym(m),Tm(m),p(m), '-r','LineWidth',2); hold off
axis tight, view(100,35)

%% Water production, with breakthrough indicated by red line
figure,
p  = cellfun(@(x) abs(x(1).qW)', wellSols,'UniformOutput',false);
p  = vertcat(p{:})*day;
j  = sum(~(p>1e-3));
m  = sub2ind(size(p),j,1:numel(W(1).cells));
surf(Ym,Tm,p); shading interp; colormap(parula(20));
hold on; plot3(Ym(m),Tm(m),p(m), '-r','LineWidth',2); hold off
axis tight, view(45,35)

%% Plot surface rates etc using GUI
plotWellSols(wellSols, cumsum(dT));

%% Capillary Equilibrium within Vertical Columns
% We watch a sharp interface form into a capillary fringe in the vertical
% direction
mrstModule add incomp

%% Grid, permeability, and fluid object
gravity reset on
exmpl = 1;
G  = computeGeometry(cartGrid([20, 1, 40], [100 1 100]));
if exmpl==1
   perm  = @(x) (350*x/100 + 50).*milli*darcy;
   rock  = makeRock(G, perm(G.cells.centroids(:,1)), .1);
   dT    = .01;
   histb = false;
else
   load rndseed.mat; rng(S);
   b = log(milli*darcy);
   a = (log(darcy)-b)/(.4 - .05);
   p = gaussianField(G.cartDims, [0.05 0.4], [3 1 11], 4.5);
   K = exp(a*(p-.05)+b);
   rock = makeRock(G, K(:), p(:));
   dT = .1;
   histb = true;
end
hT = computeTrans(G, rock);

fluid = initSimpleFluidJfunc('mu' , [0.30860, 0.056641]*centi*poise, ...
      'rho', [ 975.86,  686.54]*kilogram/meter^3, ...
      'n' , [      2,       2], ...
      'surf_tension',1*barsa/sqrt(mean(rock.poro)/(mean(rock.perm))),...
      'rock',rock);

%% Initial data
t = 0;
xr = initResSol(G, 100.0*barsa, 0.0);
xr.s(G.cells.centroids(:,3)>50) = 1.0;

%% Plot permeability and prepare for saturation
clf, set(gcf,'Position',[0   450  1280   370]);
cax1 = subplot(1,3,1);
colormap(cax1, parula);
if exmpl==2
    plotCellData(G,log10(rock.perm),'EdgeColor','none'); view(0,0), axis tight
    [h,az] = colorbarHist(log10(rock.perm),[-14 -12],'South');
    set(h,'XTick',-14:-12,'XTickLabel',{'10', '100', '1000'});
else
    K = convertTo(rock.perm,milli*darcy);
    plotCellData(G,K,'EdgeColor','none'); view(0,0), axis tight
    [h,az] = colorbarHist(K,[50 400],'South',40);
end
set(az,'Position',get(az,'Position')-[0 0 0 .02]);

cax2 = subplot(1,3,2);
colormap(cax2, [zeros(128,1) linspace(.8,0,128)' linspace(0,.7,128)']);

%% Time loop
dt = dT*[1 1 2 2 3 3 4 4 repmat(5,[1,96])]*year;
dt = [dt(1).*sort(repmat(2.^-[1:5 5],1,1)) dt(2:end)];
s  = xr.s(:,1); 
for k = 1 : numel(dt)
   xr = incompTPFA(xr, G, hT, fluid);
   xr = implicitTransport(xr, G, dt(k), rock, fluid);
   t  = t+dt(k);
     
   % Plot solution
   cla
   plotCellData(G,xr.s(:,1),'EdgeColor','none');
   view(0,0); axis tight; caxis([0 1]);
   title(sprintf('time: %.1f yrs', t/year));
   drawnow
   
   ds = norm(s - xr.s(:,1),inf);
   if ds<1e-4, break, end
   fprintf('%e\n',ds);
   s = xr.s(:,1);
   
end

%% Draw plot of pc versus S
subplot(1,3,3)
plot(xr.s, fluid.pc(xr)/barsa,'o');
%% Well-Pair Diagnostics for 3D Subset of SPE10
% In this example, we show how one can use static tracer partition to
% visualize drainage and flooded volumes and compute well-pair diagnostics
% such as volumes, well-allocation factors, etc. We also show how to
% subdivide wells into multiple segments and use this to study the flow
% patterns in more detail.

mrstModule add diagnostics spe10 incomp coarsegrid
close all

%% Set up problem
% As our example, we consider a subsample of Model 2 from the 10th SPE
% Comparative Solution Project, but with a new well pattern that has two
% injectors in the center of the reservoir.

% Grid
cartDims = [  60,  220,  20];
physDims = [1200, 2200, 2*cartDims(end)] .* ft();   % ft -> m
G  = cartGrid(cartDims, physDims);
G  = computeGeometry(G);

% Rock
rock = getSPE10rock(1:cartDims(end));
rock.poro = max(rock.poro, 1e-4);

% Wells
W = [];
wname    = {'P1', 'P2', 'P3', 'P4', 'I1', 'I2'};
wtype    = {'bhp', 'bhp', 'bhp', 'bhp', 'bhp', 'bhp'};
wtarget  = [200,   200,   200,   200,   500,   500  ] .* barsa();
wrad     = [0.125, 0.125, 0.125, 0.125, 0.125, 0.125] .* meter;
wloc     = [  1,   60,     1,   60,  20, 40;
              1,    1,   220,  220, 130, 90];
for w = 1 : numel(wtype),
   W = verticalWell(W, G, rock, wloc(1,w), wloc(2,w), 1 : cartDims(end), ...
                    'Type', wtype{w}, 'Val', wtarget(w), ...
                    'Radius', wrad(w), 'Name', wname{w}, ...
                    'InnerProduct', 'ip_tpf');
end

% Fluid: single-phase flow
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);

%% Show model setup
% We start by visualizing the reservoir and well pattern. We add a colorbar
% with a histogram and fiddle somewhat with the graphics handles to get a
% nice look of the overall figure
fig1=figure('position',[450 450 750 350]);
plotCellData(G,rock.poro, 'EdgeColor','k','EdgeAlpha',.05);
plotWell(G,W,'height',2,'FontSize',14); axis tight;
set(gca,'dataaspect',[1 1 0.06]), view(-60,15); axis off

[hc,hh]=colorbarHist(rock.poro,[.0 .45],'South', 80);
pos=get(hc,'Position'); set(hc,'Position',max(pos - [-.03 0 .25 .02], 0),'FontSize',12);
pos=get(hh,'Position'); set(hh,'Position',max(pos - [-.03 0.02 .25 0.03], 0));
set(gca,'Position',[.13 .17 .775 .785])

%% Compute flow solution and diagnostic quantities
% As usual, we only solve a single-phase flow problem to get the flow field
% needed for the diagnostics
rS = initState(G, W, 0);
T  = computeTrans(G, rock);
rS = incompTPFA(rS, G, T, fluid, 'wells', W);
D = computeTOFandTracer(rS, G, rock, 'wells', W);

%% Sweep volumes for I1 and I2
% We start by showing the sweep volumes for the two injectors. Here, we observe
% three regions and not the two we expected. The third, and very small,
% region corresponds to a section of the reservoir that is almost
% impermeable and hence will not be swept by any of the injetors.
clf(fig1);
plotCellData(G,D.ipart,'EdgeColor','k','EdgeAlpha',.05);
plotWell(G,W,'height',2,'FontSize',14); axis tight; title('Flooded volumes');
set(gca,'dataaspect',[1 1 0.06]), view(-60,15); axis off

%% Drainage volumes for P2 to P4
% Next, we show the drainage volumes for the producers. To be able to
% look inside the model, we exclude the drainage region for P1, which is
% closest to the view point.
clf(fig1);
Gbb = cartGrid([1 1 1]);
Gbb.nodes.coords = bsxfun(@mtimes,Gbb.nodes.coords,max(G.nodes.coords));
plotCellData(G,D.ppart,D.ppart>1,'EdgeColor','k','EdgeAlpha',.05);
plotGrid(Gbb,'FaceColor','none');
plotWell(G,W,'height',2,'FontSize',14); axis tight; title('Drainage volumes');
set(gca,'dataaspect',[1 1 0.06]), view(-60,15); axis off

%% Show refined partition
% When using a majority vote to determine drainage and sweep regions, we
% disregard the fact that there are regions that are influenced by more
% than one tracer. In this visualization, we will blend in a gray color to
% signify that some regions are influenced by more than one tracer
clf(fig1),
plotTracerBlend(G, D.ppart, max(D.ptracer, [], 2), 'EdgeColor','k','EdgeAlpha',.05)
plotWell(G,W,'height',2,'FontSize',14); axis tight;
set(gca,'dataaspect',[1 1 0.06]), view(-60,15); axis off

%% Well pairs
% Having established the injection and tracer partitions, we can identify
% well pairs and compute the pore volumes of the regions that can be
% associated with each well pair.
fig2=figure;
WP = computeWellPairs(rS, G, rock, W, D);
pie(WP.vols, ones(size(WP.vols)));
legend(WP.pairs,'location','Best');

figure(fig1); clf
p = compressPartition(D.ipart + D.ppart*max(D.ipart))-1;
plotCellData(G,p,p>0,'EdgeColor','k','EdgeAlpha',.05);
plotWell(G,W,'height',2,'FontSize',14); axis tight;
set(gca,'dataaspect',[1 1 0.06]), view(-60,15); axis off

%% Allocation factors for well pairs
% We show a bar plot of well allocation factors for each completion of the
% wells as a function of the depth of the completion. The allocation factor
% is defined as the normalized, cummulative flux in/out of a well from
% bottom and up.
% To investigate how the volumetric connections affect the flow in and out
% of wells, we can look at well-allocation factors, which are defined as
% the cumulative flux in/out of a well from toe to heel (here: bottom to
% top perforation). First, we compute the flux allocation manually for the
% two injectors
figure(fig1); clf
for i=1:numel(D.inj)
   subplot(1,numel(D.inj),i);
   alloc = cumsum(flipud(WP.inj(i).alloc),1);
   barh(flipud(WP.inj(i).z), alloc,'stacked'); axis tight
   lh=legend(W(D.prod).name,'Location','SouthEast');
   set(gca,'YDir','reverse');   title(W(D.inj(i)).name);
end
%%
% Then we use a library function to compute and visualize the
% well-allocation factors for all the wells in the model
figure; set(gcf,'Position',[10 70 600 760]);
plotWellAllocationPanel(D, WP);

%% Look at individual completions
% To look more closely at the performance of the different completions
% along the well path, we can divide the completion intervals into bins and
% assign a corresponding set of pseudo wells for which we recompute flow
% diagnostics. As an example, we split the completions of I1 into three
% bins and the completions of I2 into four bins.
[rSp,Wp] = expandWellCompletions(rS,W,[5, 3; 6, 4]);
Dp = computeTOFandTracer(rSp, G, rock, 'wells', Wp);

%% Injector I1
figure(fig1); clf
plotCellData(G, Dp.ipart,(Dp.ipart>0) & (Dp.ipart<4), 'EdgeColor','k','EdgeAlpha',.05);
plotWell(G,W,'height',2,'FontSize',14); axis tight;
plotGrid(Gbb,'FaceColor','none');
set(gca,'dataaspect',[1 1 0.06]), view(-60,15); axis off

%% Injector I2
figure(fig1); clf
plotCellData(G, Dp.ipart,(Dp.ipart>3) & (Dp.ipart<8), 'EdgeColor','k','EdgeAlpha',.05);
plotWell(G,W,'height',2,'FontSize',14); axis tight;
plotGrid(Gbb,'FaceColor','none');
set(gca,'dataaspect',[1 1 0.06]), view(120,15); axis off

%% 
% We end the example by computing the fraction of the sweep region for I2
% that can be attributed to the different well segments. To this end, we
% first recompute well-pair regions for all well segments and then use
% accumarray to sum all well pairs that involve segments from I2.
figure(fig2); clf
WPp = computeWellPairs(rSp, G, rock, Wp, Dp);
avols = accumarray(WPp.pairIx(:,1),WPp.vols);
h = pie(avols(4:end)); set(h(2:2:end),'FontSize',16);

%% Copyright notice

% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>

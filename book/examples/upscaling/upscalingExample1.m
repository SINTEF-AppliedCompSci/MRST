%% Illustrate use of flow diagnostics to verify upscaling
% In this example we explain how to use match in cumulative well-allocation
% factors to verify the quality of an upscaling.  In each well completion,
% the well-allocation factor is the percentage of the flux in/out of the
% completion that can be attributed to a pair of injection and production
% wells. The script can run two types of upscaling, if use_trans=false, we
% use a simple flow-based method to upscale permeability, while if
% use_trans=true, we use a more comprehensive method to upscale the
% transmissibilities and the well indices.
mrstModule add agglom upscaling coarsegrid diagnostics incomp;

use_trans = false;

%% Make model
% We make a small model that consists of two different facies with
% contrasting petrophysical properties. An injector and a producer are
% placed diagonally opposite of each other and completed mainly in the
% high-permeable part of the model.
G  = computeGeometry(cartGrid([40 20 15]));
K1 = gaussianField(G.cartDims, [200 2000]);
p1 = K1(:)*1e-4 + .2;
K1 = K1(:)*milli*darcy;
K2 = gaussianField(G.cartDims, [10 500]);
p2 = K2(:)*1e-4 + .2;
K2 = K2(:)*milli*darcy;

rad1 = G.cells.centroids(:,1).^2 + .5*G.cells.centroids(:,2).^2 ...
   + (G.cells.centroids(:,3)-2).^2;
rad2 = .5*(G.cells.centroids(:,1)-40).^2 + 2*G.cells.centroids(:,2).^2 ...
   + 2*(G.cells.centroids(:,3)-2).^2;

ind = ((rad1>600) & (rad1<1500)) | ((rad2>700) & (rad2<1400));
rock.perm = K2(:);
rock.perm(ind) = K1(ind);
rock.perm = bsxfun(@times, rock.perm, [1 1 1]);
rock.poro = p2(:);
rock.poro(ind) = p1(ind);
pv = poreVolume(G, rock);

W = verticalWell([], G, rock, 4, 17, 4:15, 'Comp_i', [1 0], ...
   'Type', 'rate', 'Val', 0.2*sum(pv)/year, 'Name', 'I');
W = verticalWell(W,  G, rock, 35, 3, 1:10, 'Comp_i', [0 1], ...
   'Type', 'rate', 'Val', -0.2*sum(pv)/year, 'Name', 'P');

figure(1); clf,
set(gcf,'Position', [860 450 840 400],'PaperPositionMode','auto');
subplot(1,2,1);
plotCellData(G,rock.poro,'EdgeAlpha',.5); view(3);
plotWell(G,W,'Color','k'); axis off tight
pax = [min(rock.poro) max(rock.poro)];
[hc,hh] = colorbarHist(rock.poro, pax,'West');
pos = get(hh,'Position'); set(hh,'Position', pos - [.01 0 0 0]);
pos = get(hc,'Position'); set(hc,'Position', pos - [.03 0 0 0]);
set(gca,'Position',[.12 .075 .315 .9])

%% Compute flow on the fine-scale model
% Transmissibility
hT = computeTrans(G, rock);
trans = 1 ./ accumarray(G.cells.faces(:,1), 1 ./ hT, [G.faces.num, 1]);

% Fluid model
dfluid = initSimpleFluid('mu' , [   1,  10]*centi*poise     , ...
                         'rho', [1014, 859]*kilogram/meter^3, ...
                         'n', [2 2]);
gravity reset off;

% Pressure solution
xd  = initState(G, W, 100*barsa);
nw = numel(W);
xd  = incompTPFA(xd, G, trans, dfluid, 'wells', W, 'use_trans', true);


%% Flow diagnostics on fine-scale model
% To better reveal the communication in the reservoir, we subdivide each
% well into two segments, an upper and a lower segments, so that we
% altogether will have four well-pairs whose well-allocation factors can be
% used to verify the quality of the upscaling
nsegment=2;
[xd,Wdf] = expandWellCompletions(xd, W,[(1:nw)' repmat(nsegment,nw,1)]);
Df  = computeTOFandTracer(xd, G, rock, 'wells', Wdf);
WPf = computeWellPairs(xd, G, rock, Wdf, Df);

%% Coarse-scale solution
% We make a 5x5x15 coarse grid, in which we have chosen to keep the
% layering of the fine-scale model to simplify the comparison of
% well-allocation factors. Transmissibilities and well indices are upscaled
% using two slightly different methods: for the transmissibilities we use
% global generic boundary conditions and on each coarse face use the
% solution that has the largest flux orthogonal to the face to compute the
% upscaled transmissibility. For the wells, we use specific well
% conditions and use least squares for the flux.
p  = partitionCartGrid(G.cartDims, [5 5 15]);
CG = coarsenGeometry(generateCoarseGrid(G, p));
crock = convertRock2Coarse(G, CG, rock);
if use_trans  %#ok<*UNRCH>
   [~,CTrans] = upscaleTrans(CG, hT, 'match_method', 'max_flux', ...
                          'bc_method', 'bc_simple');
   [~,~,WC]   = upscaleTrans(CG, hT, 'match_method', 'lsq_flux', ...
                          'bc_method', 'wells_simple', 'wells', W);
else
   crock.perm = upscalePerm(G, CG, rock);
end
figure(1); subplot(1,2,2); cla
plotCellData(CG,crock.poro,'EdgeColor','none');
plotFaces(CG, boundaryFaces(CG),...
   'EdgeColor', [0.4 0.4 0.4],'EdgeAlpha',.5, 'FaceColor', 'none'); view(3);
plotWell(G,W,'Color','k'); axis off tight
[hc,hh]=colorbarHist(crock.poro, pax,'East');
pos = get(hh,'Position'); set(hh,'Position', pos + [.01 0 0 0]);
pos = get(hc,'Position'); set(hc,'Position', pos + [.02 0 0 0]);
set(gca,'Position',[.56 .075 .315 .9])

%%
% If permeability upscaling: to assess the effect of permeability upscaling
% only, we project the upscaled permeability back to the fine grid
if ~use_trans
   crock.poro = crock.poro(p);
   crock.perm = crock.perm(p,:);
   CG = G;
   WC = W;
   CTrans = computeTrans(CG, crock);
   CTrans = 1 ./ accumarray(CG.cells.faces(:,1), 1 ./ CTrans, [CG.faces.num, 1]);
   p = (1:G.cells.num)';
end
xd    = initState(CG, WC, accumarray(p,100*barsa)./accumarray(p,1));
xd    = incompTPFA(xd, CG, CTrans, dfluid, 'wells', WC, 'use_trans', true);

%% Flow diagnostics on the upscaled model
% Having obtained fluxes on the coarse model, we can expand the wells into
% two segments that exactly match the subdivision of in the fine-scale
% model and compute flow diagnostics.
nw    = numel(WC);
[xdc,Wdc] = expandCoarseWellCompletions(xd, WC, Wdf, p);
Dc    = computeTOFandTracer(xdc, CG, crock, 'wells', Wdc);
WPc   = computeWellPairs(xdc, CG, crock, Wdc, Dc);
nit   = numel(Df.inj);
npt   = numel(Df.prod);
nseg  = nit + npt;

%% Display bar charts of flow-allocation factors
% For each well segment, we plot bar chart of the cumulative flux in/out of
% the completions that make up the segment, from bottom to top. In the
% plots, each well segment is assigned a unique color and each bar is
% subdivided into the fraction of the total in/outflux that blongs to the
% different well-pairs the segment is part of.
%
% The plots show the allocation factors for the upper half of the injector
% (I:1, upper-left plot), the lower part of the injector (I:2, upper-right
% plot), the upper part of the producer (P:1, lower-left plot), and the
% lower part of the producer (P:2, lower-right plot). Looking at the bar
% chart for I:1, we see that the majority of the flux from this injector
% is colored yellow and hence goes to P:1, which corresponds to the upper
% part of the producer.
figure(2); set(gcf,'Position', [860 450 840 400]); clf
plotWellAllocationComparison(Dc, WPc, Df, WPf);
dy = [-.05 -.05 -.025 -.025];
dx = [-.5 0 -.5 0 0];
for i=1:4
   subplot(2,2,i);
   pos = get(gca,'Position');
   pos = pos + [-.05 dy(i) .075 .075];
   set(gca,'Position',pos,'XTick',[],'YTick',[]);
end
cmap = jet(nseg); cmap = 0.6*cmap + .4*ones(size(cmap)); colormap(cmap);


%% Display the partition
% Last we show the cells that belong to the four different well segments
% and show the corresponding injector/producer partitions.
%
% We can now compare the partitions with the well-allocation factors in
% Figure 2. Let us take the allocation for I:2 as an example. In the
% upper-right plot of Figure 2 we see that almost all the flux from this
% injector goes to P:2. Looking at the cyan region in the left plot of
% Figure 1 confirms this: this region is hardly in contact with the yellow
% segment of the producer. The blue region, on the other hand, contains
% both the upper segment of the producer (P:1, yellow color) and parts of
% the lower segment (P:2, red color). Likewise, the red volume in the
% right-hand plot of Figure 3 covers both the lower segment of the injector
% (I:2, cyan) and parts of the upper segment (I:1, blue). In the bars of
% P:2 to the lower-right in Figure 2, we therefore see that relatively
% large portion of the bars are in blue color.
oG  = computeGeometry(cartGrid([1 1 1],[40 20 15]));
figure(3); clf
set(gcf,'Position', [860 450 840 310],'PaperPositionMode','auto');
subplot(1,2,1);
for i=1:nit
   plotCellData(G,Df.ipart,Df.ipart==i,'FaceAlpha',.3,'EdgeAlpha',.2);
end
plotGrid(oG,'FaceColor','none');
for i=1:nseg,
   plotGrid(G,Wdf(i).cells,'FaceColor',cmap(i,:));
end
plotWell(G,W,'Color','k');
view(-40,10); axis tight off; caxis([1 nit+npt]);

subplot(1,2,2);
for i=1:npt
   plotCellData(G,Df.ppart+nit,Df.ppart==i,'FaceAlpha',.4,'EdgeAlpha',.2);
end
plotGrid(oG,'FaceColor','none');
for i=1:nseg,
   plotGrid(G,Wdf(i).cells,'FaceColor',cmap(i,:));
end
plotWell(G,W,'Color','k');
view(-40,10); axis tight off; caxis([1 nit+npt]);
colormap(cmap);

%% Arithmetic-harmonic averaging
% To illustrate the implementation of this type of averaging, we consider a
% single layer of the SPE10 model
fine = [45 60];
G    = cartGrid(fine, fine);
G    = computeGeometry(G);
rock = SPE10_rock(1:fine(1),1:fine(2),46);
rock.perm = convertFrom(rock.perm(:,1:2), milli*darcy);

% The coarse grid
coarse = [15 15];
q = partitionUI(G,coarse);
CG = generateCoarseGrid(G, q);

%%
%
vol = G.cells.volumes;
for i=1:size(rock.perm,2)
   crock1.perm(:,i) = accumarray(q,vol.*rock.perm(:,i))./accumarray(q,vol);
   
   crock2.perm(:,i) = accumarray(q,vol)./accumarray(q,vol./rock.perm(:,i));

   dims = fine; dims(i)=coarse(i);
   qq = partitionUI(G, dims);
   K = accumarray(qq,vol)./accumarray(qq,vol./rock.perm(:,i));
   crock3.perm(:,i) = accumarray(q,K(qq).*vol)./accumarray(q,vol);
end

%%
px = log10([min(rock.perm(:,1)) max(rock.perm(:,1))]);
clf
subplot(1,4,1); 
plotCellData(G,log10(rock.perm(:,2)),'EdgeColor','none');
plotGrid(CG,'FaceColor','none','EdgeColor','k'); caxis(px); axis tight off

subplot(1,4,2);
plotCellData(G,log10(crock1.perm(q,2)),'EdgeColor','none');
plotGrid(CG,'FaceColor','none','EdgeColor','k'); caxis(px); axis tight off

subplot(1,4,3);
plotCellData(G,log10(crock2.perm(q,2)),'EdgeColor','none');
plotGrid(CG,'FaceColor','none','EdgeColor','k'); caxis(px); axis tight off

subplot(1,4,4);
plotCellData(G,log10(crock3.perm(q,2)),'EdgeColor','none');
plotGrid(CG,'FaceColor','none','EdgeColor','k'); caxis(px); axis tight off
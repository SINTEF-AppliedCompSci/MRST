%% Averaging of SPE10 layer
% To illustrate the implementation of various types of averaging
% techniques, we consider a single layer of the SPE10 model and compare
% the effective permeabilities obtained by arithmetic, harmonic, and
% harmonic-arithmetic averaging.
mrstModule add spe10 coarsegrid;

%% Setup model
fine = [40 60];
G    = cartGrid(fine, fine);
G    = computeGeometry(G);
rock = getSPE10rock(1:fine(1),1:fine(2),46);
rock.perm = rock.perm(:,1:2);

% The coarse grid
coarse = [10 15];
q = partitionUI(G,coarse);
cG = cartGrid(coarse,fine);

%% Compute the averages
% The arithmetic and harmonic averages are straightforward. To compute the
% harmonic-arithmetic average, we  temporary partition that coincides with
% the coarse grid in along the given axial direction and with the original
% fine grid in the other directions. Then we map the averaged values back
% onto the fine grid and perform a standard arithmetic averaging.
clear crock*
vol = G.cells.volumes;
for i=1:size(rock.perm,2)
   crock1.perm(:,i) = accumarray(q,vol.*rock.perm(:,i))./accumarray(q,vol);
   
   crock2.perm(:,i) = accumarray(q,vol)./accumarray(q,vol./rock.perm(:,i));

   dims = G.cartDims; dims(i)=coarse(i);
   qq = partitionUI(G, dims);
   K = accumarray(qq,vol)./accumarray(qq,vol./rock.perm(:,i));
   crock3.perm(:,i) = accumarray(q,K(qq).*vol)./accumarray(q,vol);
end

%% Plot the resulting effective permeabilities
% To improve the visualization we enhance the colorbar by adding a small
% histogram of permeability values
px = log10([min(rock.perm(:,1)) max(rock.perm(:,1))]);
clf, set(gcf,'Position',[500 400 880 410]);
subplot(1,4,1); 
plotCellData(G,log10(rock.perm(:,2)),'EdgeColor','none');
plotGrid(cG,'FaceColor','none','EdgeColor','k'); caxis(px); axis tight off
title('Fine scale');

subplot(1,4,2);
plotCellData(cG,log10(crock1.perm(:,2)),'EdgeColor','k');
caxis(px); axis tight off
title('Arithmetic');

subplot(1,4,3);
plotCellData(cG,log10(crock2.perm(:,2)),'EdgeColor','k');
caxis(px); axis tight off
title('Harmonic');

subplot(1,4,4);
plotCellData(cG,log10(crock3.perm(:,2)),'EdgeColor','k');
caxis(px); axis tight off
title('Harmonic-arithmetic');

perm{1} = log10(rock.perm(:,2));
perm{2} = log10(crock1.perm(:,2));
perm{3} = log10(crock2.perm(:,2));
perm{4} = log10(crock3.perm(:,2));

for i=1:4
   subplot(1,4,i);
   [h,ax] = colorbarHist(perm{i},[-18 -10],'South');
   set(h,'XTick',-16:2:-12,'XTickLabel',{'.1', '10', '1000'});
end
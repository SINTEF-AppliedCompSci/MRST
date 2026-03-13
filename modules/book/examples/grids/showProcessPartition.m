load(fullfile('data','showProcessPartition.mat'));

%% Plot grid
figure, plotGrid(G,c,'FaceColor',0.7*col+[.3 .3 .3]);
text(G.cells.centroids(c,1), G.cells.centroids(c,2),...
    num2str((1:numel(c))'),'HorizontalAlignment','center','FontSize',20);
axis off tight

%% Plot original adjacency matrix
figure,
val = ones(size(c));
val(p(1:r(2)-1))=1; val(p(r(2):r(3)-1))=2;
hold on;
i1 = val(ii)==1; plot(ii(i1),jj(i1),'or','MarkerSize',6,'MarkerFaceColor','r');
i2 = val(ii)==2; plot(ii(i2),jj(i2),'ob','MarkerSize',6,'MarkerFaceColor','b');
set(gca,'YDir','reverse','FontSize',14); axis square tight
hold off

%% Plot permuted adjacency matrix
figure,
n(p) = 1:numel(c);
hold on;
plot(n(ii(i1)),n(jj(i1)),'or','MarkerSize',6,'MarkerFaceColor','r');
plot(n(ii(i2)),n(jj(i2)),'ob','MarkerSize',6,'MarkerFaceColor','b');
plot([12.5 12.5 NaN .5 24.5],[.5 24.5 NaN 12.5 12.5],'--k');
set(gca,'YDir','reverse','FontSize',14); axis square tight
hold off

clc; clear all; close all;

xMax = 2; yMax = 2;
G = cartGrid([1,1], [xMax, yMax]);
G.nodes.coords = bsxfun(@minus, G.nodes.coords,[1,1]);
epsilon = .6;
G.nodes.coords(G.nodes.coords(:,1) == 1 & G.nodes.coords(:,2) == 1,:) ...
    = [1,1+epsilon];

G = computeGeometry(G);
plotGrid(G, 'facealpha', .2)
set(gcf, 'defaulttextinterpreter', 'latex')

axis equal
axis([-1.5, 1.5, -1.5, 1+epsilon+.5])

set(gca,'XTick',[-1 1] );
set(gca, 'XTickLabels', {'',''})
set(gca,'YTick',[-1 1, 1+epsilon] );
set(gca, 'YTickLabels', {'','',''})
my_xticklabels([-1,1], {'-1','1'},0,-.1)
my_yticklabels([-1,1,1+epsilon], {'-1','-1','1+$\varepsilon$'},.1,-.2)

xlabel('$x$'); ylabel('$y$');

cut = 0;
ps = get(gcf, 'Position');
ratio = 1;
paperWidth = 10;
paperHeight = paperWidth*ratio;
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'papersize', [paperWidth-cut paperHeight-cut]);
set(gcf, 'PaperPosition', [-cut    -cut   paperWidth+cut paperHeight+cut]);

print(gcf, '-dpdf', '../../tex/thesis/fig/Keps.pdf');
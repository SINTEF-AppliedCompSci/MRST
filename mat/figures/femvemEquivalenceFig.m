clc; clear all; close all;
addpath('../');

add = .5;
lineWidth = 1.5;
faceAlpha = .2;
mrkSzBig = 20;
mrkSzSml = 10;
cut = 4;
addTxt = .25;
azel = [-30,-15];

%%

P = [-1,-1; 1,-1; -1,1];
A = [2,2;-3,1];

P = P*A;

patch(P(:,1), P(:,2), 'y', 'facealpha', faceAlpha)
axis equal off;

cut =1;
ps = get(gcf, 'Position');
ratio = 1;
paperWidth = 10;
paperHeight = paperWidth*ratio;
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'papersize', [paperWidth-cut paperHeight-cut]);
set(gcf, 'PaperPosition', [-cut    -cut   paperWidth+cut paperHeight+cut]);

print(gcf, '-dpdf', '../../tex/thesis/fig/tri.pdf');

%%
clf;

P = [-1, -1, -1; 1 -1 -1; -1 1 -1;1 1 1];
A = [1,-2,3; 3,4,-1; 1,2,1];
P = P*A;

ii = { 1:3, ...
       2:4, ...
       [1,2,4], ...
       [1,2,4]};

for i = 1:numel(ii)
    patch(P(ii{i},1), P(ii{i},2), P(ii{i},3), 'y', 'facealpha', .2);
end
axis equal off
view([56,-21])

%%

cut =1;
ps = get(gcf, 'Position');
ratio = 1;
paperWidth = 10;
paperHeight = paperWidth*ratio;
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'papersize', [paperWidth-cut paperHeight-cut]);
set(gcf, 'PaperPosition', [-cut    -cut   paperWidth+cut paperHeight+cut]);

print(gcf, '-dpdf', '../../tex/thesis/fig/tet.pdf');


%%
clf;

G = cartGrid([1,1]);

A = [-1,2;-3,1];
G.nodes.coords = G.nodes.coords*A;

plotGrid(G, 'facealpha', faceAlpha);

axis equal off

%%

ps = get(gcf, 'Position');
ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
paperWidth = 10;
paperHeight = paperWidth*ratio - cut;
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'papersize', [paperWidth paperHeight]);
set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);

print(gcf, '-dpdf', '../../tex/thesis/fig/quad.pdf');

%%
clf;

G = cartGrid([1,1,1]);

A = [3,2,1; 1,2,6; 1,2,1];
G.nodes.coords = G.nodes.coords*A;

plotGrid(G, 'facealpha', faceAlpha);
view([18,38])

% xMax = max(G.nodes.coords(:,1)); xMin = min(G.nodes.coords(:,1));
% yMax = max(G.nodes.coords(:,2)); yMin = min(G.nodes.coords(:,2));
% zMax = max(G.nodes.coords(:,3)); zMin = min(G.nodes.coords(:,3));
% axis([xMin, xMax, yMin, yMax, zMin, zMax]);
axis equal off


%%

cut =2;
ps = get(gcf, 'Position');
ratio = 1;
paperWidth = 10;
paperHeight = paperWidth*ratio;
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'papersize', [paperWidth-cut paperHeight-cut]);
set(gcf, 'PaperPosition', [-cut    -cut   paperWidth+cut paperHeight+cut]);

print(gcf, '-dpdf', '../../tex/thesis/fig/pol.pdf');

%%
clf;

P = [-1,-1,-1;
      1,-1,-1;
     -1, 1,-1;
     -1,-1, 1;
      1,-1, 1;
     -1, 1, 1];
 
 A = [3,2,3; 1,2,6; 1,2,1];
 P = P*A;

 ii = {[1,2,5,4], ...
       [2,3,6,5], ...
       [1,3,6,4], ...
       [1:3]    , ...
       [4:6]   };

 for i = 1:numel(ii)
    patch(P(ii{i},1), P(ii{i},2), P(ii{i},3), 'y', 'facealpha', .2);
 end
 
 axis equal off
 view([2,-65])
 
 %%
cut =3;
ps = get(gcf, 'Position');
ratio = 1;
paperWidth = 10;
paperHeight = paperWidth*ratio;
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'papersize', [paperWidth-cut paperHeight-cut]);
set(gcf, 'PaperPosition', [-cut    -cut   paperWidth+cut paperHeight+cut]);

print(gcf, '-dpdf', '../../tex/thesis/fig/pri.pdf');

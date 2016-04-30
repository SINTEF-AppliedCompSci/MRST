clc; clear; close all;

addpath('../')
addpath('../VEM2D/')

ntnuBlue = [0,61,242]/255;

G = cartGrid([2,2]);
k = 1;
f = @(X) X(:,1);
addy = .5;
lineWidth = 1.5;
faceAlpha = .2;
mrkSzBig = 20;
mrkSzSml = 10;
cut = 4;
addTxt = .25;

X = G.nodes.coords;
X = bsxfun(@minus, X,[1,1]);
G.nodes.coords = X;

alpha = 1;

G = computeVEM2DGeometry(G);

edgeNodes = G.faces.nodes;
nN = size(edgeNodes,1)/2;
Xe = G.nodes.coords(edgeNodes,:);
Xe = reshape(Xe,2,[]);


fig1 = figure();
set(fig1,'defaulttextinterpreter','latex');
plotGrid(G, 'FaceAlpha', faceAlpha);
hold on
plot(Xe(:,1:nN), Xe(:,nN+1:2*nN),'k')
line([-(1+addy) 1+addy], [0,0], 'LineWidth', lineWidth);
line([0,0], [-(1+addy) 1+addy], 'LineWidth', lineWidth);
X1 = X((X(:,1)).^2 + (X(:,2)).^2 <= 1,:);
fig = gcf;
plot(X1(:,1), X1(:,2), 'ok', 'MarkerFaceColor', 'r')
text(1+.3, .1,'$x$', 'interpreter', 'latex')
text(.1, 1+.3, '$y$', 'interpreter', 'latex')
text(-.5,-1.2,'$h_x$', 'interpreter', 'latex')
text(-1.3,-.5,'$h_y$', 'interpreter', 'latex')

set(gca,'XTick',[-1 0 1] );
set(gca,'XTickLabel',{'', '', ''})
set(gca,'YTick',[-1 0 1] );
set(gca,'YTickLabel',{'', '', ''} )
yWidth = .05;
xWidth = .35;
addy = .1;
addx = .2;
text(-1-xWidth/2       , -1.5-yWidth-addy, '$i-1$');
text(0-xWidth/6        , -1.5-yWidth-addy, '$i$'  );
text(1-xWidth/2        , -1.5-yWidth-addy, '$i+1$');
text(-1.5-xWidth-addx  , -1              , '$j-1$');
text(-1.5-xWidth/3-addx,  0              , '$j$'  );
text(-1.5-xWidth -addx ,  1              , '$j+1$');

axis equal
axis([-1.5, 1.5, -1.5 1.5])

%%

ps = get(gcf, 'Position');
ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
paperWidth = 10;
paperHeight = paperWidth*ratio - cut;
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'papersize', [paperWidth paperHeight]);
set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);

print(gcf, '-dpdf', '../../tex/thesis/fig/2Dstencil1.pdf');

% axis equal


fig2 = figure();
set(fig2,'defaulttextinterpreter','latex');
plotGrid(G, 'FaceAlpha', faceAlpha);
hold on
plot(Xe(:,1:nN), Xe(:,nN+1:2*nN),'k')
addy = .5;
line([-(1+addy) 1+addy], [-(1+addy),1+addy], 'LineWidth', lineWidth);
line([1+addy,-(1+addy)], [-(1+addy) 1+addy], 'LineWidth', lineWidth);
X1 = [X((X(:,1)).^2 + (X(:,2)).^2 == 2,:); [0,0]];
% plot(X(:,1), X(:,2), '.k', 'MarkerSize', mrkSzSml);
plot(X(:,1), X(:,2), 'ok', 'MarkerFaceColor', 'r')
text(1-.1, 1+.3,'$d_1$')
text(-(1-.1), 1+.3, '$d_2$')

set(gca,'XTick',[-1 0 1] );
set(gca,'XTickLabel',{'', '', ''})
set(gca,'YTick',[-1 0 1] );
set(gca,'YTickLabel',{'', '', ''} )
yWidth = .05;
xWidth = .35;
addy = .1;
addx = .2;
text(-1-xWidth/2       , -1.5-yWidth-addy, '$i-1$');
text(0-xWidth/6        , -1.5-yWidth-addy, '$i$'  );
text(1-xWidth/2        , -1.5-yWidth-addy, '$i+1$');
text(-1.5-xWidth-addx  , -1              , '$j-1$');
text(-1.5-xWidth/3-addx,  0              , '$j$'  );
text(-1.5-xWidth -addx ,  1              , '$j+1$');


axis equal
axis([-1.5, 1.5, -1.5 1.5])

ps = get(gcf, 'Position');
ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
paperWidth = 10;
paperHeight = paperWidth*ratio - cut;
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'papersize', [paperWidth paperHeight]);
set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);

print(gcf, '-dpdf', '../../tex/thesis/fig/2Dstencil2.pdf');
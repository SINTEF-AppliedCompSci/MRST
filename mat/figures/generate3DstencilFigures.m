clc; clear all; close all;

addpath('../VEM3D/')

G = cartGrid([2,2,2]);
k = 1;
f = @(X) X(:,1);
add = 1;


X = G.nodes.coords;
X = bsxfun(@minus, X,[1,1,1]);
G.nodes.coords = X;

G = computeVEM3DGeometry(G);

ntnuBlue = [0,61,242]/255;

lineWidth = 1.5;
faceAlpha = .2;
mrkSzBig = 20;
mrkSzSml = 10;
cut = 4;
azel = [-30,-15];


edgeNodes = G.edges.nodes;
nN = size(edgeNodes,1)/2;
Xe = G.nodes.coords(edgeNodes,:);
Xe = reshape(Xe,2,[]);


%%

fig1 = figure();
set(gcf, 'defaultTextInterpreter', 'latex')
plotGrid(G, 'FaceAlpha', .2);
edgeNodes = G.edges.nodes;
nN = size(edgeNodes,1)/2;
hold on
plot3(Xe(:,1:nN), Xe(:,nN+1:2*nN), Xe(:,2*nN+1:3*nN),'k')
X1 = X((X(:,1)).^2 + (X(:,2)).^2 + (X(:,3)).^2 <= 1,:);
line([-(1+add) 1+add], [0,0], [0,0], 'LineWidth', lineWidth);
line([0,0], [-(1+add) 1+add], [0,0], 'LineWidth', lineWidth);
line([0,0], [0,0], [-(1+add) 1+add], 'LineWidth', lineWidth);
plot3(X1(:,1), X1(:,2), X1(:,3), 'ok', 'MarkerFaceColor', 'r')

set(gca,'XTick',[-1 0 1] );
set(gca,'XTickLabel',{'i-1', 'i', 'i+1'} )
set(gca,'YTick',[-1 0 1] );
set(gca,'YTickLabel',{'j-1', 'j', 'j+1'} )
set(gca,'ZTick',[-1 0 1] );
set(gca,'ZTickLabel',{'k-1', 'k', 'k+1'} )
set(gca, 'FontName', 'Times-Roman', 'FontAngle', 'Oblique')

% h = my_xticklabels([-1,0,1], {'$i-1$', '$i$', '$i+1$'});
% h = my_yticklabels([-1,0,1], {'$j-1$', '$j$', '$j+1$'});

text(1+.4,0,.2,'$x$')
text(0,1+.4,.2,'$y$')
text(.15,0,1+.5,'$z$')
text(-.5,1.2,1.1,'$h_x$')
text(-1.2,.4,.8,'$h_y$')
text(-1.2,1.2,.5,'$h_z$')

axis equal
axis([-1.5 1.5 -1.5 1.5 -1.5 1.5])
view(azel)

%%

ps = get(gcf, 'Position');
ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
paperWidth = 10;
paperHeight = paperWidth*ratio - cut;
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'papersize', [paperWidth paperHeight]);
set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);

print(gcf, '-dpdf', '../../tex/thesis/fig/3Dstencil1.pdf');

%%

fig2 = figure();
set(gcf, 'defaultTextInterpreter', 'latex')
plotGrid(G,'FaceAlpha', .2);
edgeNodes = G.edges.nodes;
nN = size(edgeNodes,1)/2;
Xe = G.nodes.coords(edgeNodes,:);
Xe = reshape(Xe,2,[]);
hold on
plot3(Xe(:,1:nN), Xe(:,nN+1:2*nN), Xe(:,2*nN+1:3*nN),'k')
X2 = [X( ((X(:,2).^2 + X(:,3).^2 <= 2) & (X(:,1) == 0)),:); [-1,0,0]; [1,0,0]];

line([-(1+add) 1+add], [0,0], [0,0], 'LineWidth', lineWidth);
line([0,0],[-sqrt(1+add) sqrt(1+add)], [-sqrt(1+add) sqrt(1+add)], 'LineWidth', lineWidth);
line([0,0],[-sqrt(1+add) sqrt(1+add)], [sqrt(1+add) -sqrt(1+add)], 'LineWidth', lineWidth);
plot3(X2(:,1), X2(:,2), X2(:,3), 'ok', 'MarkerFaceColor', 'r')

set(gca,'XTick',[-1 0 1] ); %This are going to be the only values affected.
set(gca,'XTickLabel',{'i-1', 'i', 'i+1'} )
set(gca,'YTick',[-1 0 1] ); %This are going to be the only values affected.
set(gca,'YTickLabel',{'j-1', 'j', 'j+1'} )
set(gca,'ZTick',[-1 0 1] ); %This are going to be the only values affected.
set(gca,'ZTickLabel',{'k-1', 'k', 'k+1'} )
set(gca, 'FontName', 'Times-Roman', 'FontAngle', 'Oblique')
axis equal

text(1+.4,0,.2,'$x$')
text(0,1+.2,1+.3,'$\eta_x$')
text(-.05,-(1+.4),1+.25,'$\zeta_x$')

axis equal
axis([-1.5 1.5 -1.5 1.5 -1.5 1.5])
view(azel)

%%

ps = get(gcf, 'Position');
ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
paperWidth = 10;
paperHeight = paperWidth*ratio - cut;
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'papersize', [paperWidth paperHeight]);
set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);

print(gcf, '-dpdf', '../../tex/thesis/fig/3Dstencil2.pdf');

%%

fig3 = figure();
set(gcf, 'defaultTextInterpreter', 'latex')
plotGrid(G,'FaceAlpha', .2);
edgeNodes = G.edges.nodes;
nN = size(edgeNodes,1)/2;
Xe = G.nodes.coords(edgeNodes,:);
Xe = reshape(Xe,2,[]);
hold on
plot3(Xe(:,1:nN), Xe(:,nN+1:2*nN), Xe(:,2*nN+1:3*nN),'k')
% X3 = [X( X(:,1).^2 + X(:,2).^2 + X(:,3).^2 == 3, :); [0,0,0]];
line([-sqrt(1+add) sqrt(1+add)], [-sqrt(1+add) sqrt(1+add)], [-sqrt(1+add) sqrt(1+add)], 'LineWidth', lineWidth);
line([sqrt(1+add) -sqrt(1+add)], [-sqrt(1+add) sqrt(1+add)], [-sqrt(1+add) sqrt(1+add)], 'LineWidth', lineWidth);
line([-sqrt(1+add) sqrt(1+add)], [sqrt(1+add) -sqrt(1+add)], [-sqrt(1+add) sqrt(1+add)], 'LineWidth', lineWidth);
line([-sqrt(1+add) sqrt(1+add)], [-sqrt(1+add) sqrt(1+add)], [sqrt(1+add) -sqrt(1+add)], 'LineWidth', lineWidth);
plot3(X(:,1), X(:,2), X(:,3), 'ok', 'MarkerFaceColor', 'r')


set(gca,'XTick',[-1 0 1] ); %This are going to be the only values affected.
set(gca,'XTickLabel',{'i-1', 'i', 'i+1'} )
set(gca,'YTick',[-1 0 1] ); %This are going to be the only values affected.
set(gca,'YTickLabel',{'j-1', 'j', 'j+1'} )
set(gca,'ZTick',[-1 0 1] ); %This are going to be the only values affected.
set(gca,'ZTickLabel',{'k-1', 'k', 'k+1'} )
set(gca, 'FontName', 'Times-Roman', 'FontAngle', 'Oblique')
axis equal

text(1+.4,1+.2,1+.2,'$d_1$')
text(1+.4,-(1-.5),1+.2,'$d_2$')
text((.7),(.8),-(1+.1),'$d_3$')
text((1+.5),-(1-.3),-(1+.2),'$d_4$')

axis([-1.5 1.5 -1.5 1.5 -1.5 1.5])
view(azel)

%%

ps = get(gcf, 'Position');
ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
paperWidth = 10;
paperHeight = paperWidth*ratio - cut;
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'papersize', [paperWidth paperHeight]);
set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);


print(gcf, '-dpdf', '../../tex/thesis/fig/3Dstencil3.pdf');
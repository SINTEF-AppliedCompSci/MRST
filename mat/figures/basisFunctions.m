clc; clear all; close all;

addpath('../VEM2D/')
addpath('../')

load('basisGrid');

G = sortEdges(G);
G = computeVEM2DGeometry(G);
cut = 3;

%%

bEdg = find(any(G.faces.neighbors == 0,2));

n = size(P,1);
for i = 1:n

    tol = 1e-6;
    dofEdg1 = abs(G.faces.centroids(bEdg,1+l(i))          - ...
                 (G.faces.centroids(bEdg,2-l(i))*a(i) + b(i))) < tol;
    dofEdg2 = abs(G.faces.centroids(bEdg,1+l(mod(i,n)+1)) - ...
                 (G.faces.centroids(bEdg,2-l(mod(i,n)+1))*a(mod(i,n)+1) + b(mod(i,n)+1))) < tol;

    A = [P(mod(i,n)+1,:)-P(i,:);P(mod(i+1,n)+1,:)-P(i,:)]\[1;0];
    B = -P(i,:)*A;

    gD = @(X) X*A + B;

    bc = VEM2D_addBC([], G, bEdg(dofEdg1), 'pressure', gD);
    bc = VEM2D_addBC(bc, G, bEdg(dofEdg2), 'pressure', gD);
    bc = VEM2D_addBC(bc, G, bEdg(~dofEdg2 & ~dofEdg1), 'pressure', 0);

    sol = VEM2D(G, 0, bc, 2);

    figure;
    
    plotVEM2D(G, sol, 1, 'edgecolor', 'none');
    hold on
    for j = 1:6
        plot3(P(:,1), P(:,2), gD(P), 'k')
    end
    view(-37.5, 40)
    xlabel('$x$', 'interpreter', 'latex')
    ylabel('$y$', 'interpreter', 'latex')
    zlabel('$z$', 'interpreter', 'latex')
    axis equal;
    axis([0 1 0 1 0 1])
    
    set(gca,'XTick',[0 1] );
    set(gca,'YTick',[0 1] );
    set(gca,'ZTick',[0 1] );
    
    ps = get(gcf, 'Position');
    ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
    paperWidth = 10;
    paperHeight = paperWidth*ratio - cut;
    set(gcf, 'paperunits', 'centimeters');
    set(gcf, 'papersize', [paperWidth paperHeight]);
    set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);
%     dest = strcat('../../tex/thesis/fig/basis2D/Phi_', num2str(i), '.pdf');
%     print(gcf, '-dpdf', dest);
%     
    dest = strcat('../../tex/thesis/fig/basis2D/Phi_', num2str(i), '.png');
    print(gcf, '-dpng', dest, '-r600');
    
end

%%
figure;
fill(P(:,1), P(:,2), 'y', 'facealpha', .2)
hold on
xK = mean(polygonInt(G, 1:G.cells.num, @(X) X(:,1), 7)./G.cells.volumes);
yK = mean(polygonInt(G, 1:G.cells.num, @(X) X(:,2), 7)./G.cells.volumes);
plot(xK, yK, 'ok', 'MarkerFaceColor', 'r', 'markersize', 4);
text(xK*.9, yK*.9,'$\textbf{x}_K$', 'interpreter', 'latex')
plot([P(2,1),P(5,1)],[P(2,2),P(5,2)],'k--')
text(.65, .25,'$h_K$', 'interpreter', 'latex')
addFac = .1;
% for i = 1:n
%     cVec = P(i,:) - [xK, yK];
%     x = P(i,1) + cVec(1)*addFac;
%     y = P(i,2) + cVec(2)*addFac;
%     text(x,y,num2str(i),'interpreter', 'latex', 'fontsize', 8);
% end
axis equal off

%%

w = 5;
h = 5;
ps = get(gcf, 'Position');
ratio = 1;
paperWidth = 10;
paperHeight = paperWidth*ratio;
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'papersize', [paperWidth-w paperHeight-h]);
set(gcf, 'PaperPosition', [-w    -h   paperWidth+w paperHeight+h]);
% print(gcf, '-dpdf', '../../tex/thesis/fig/basis2D/BasisElement.pdf');
print(gcf, '-dpng', '../../tex/thesis/fig/basis2D/BasisElement.png', '-r600');
clc; clear all; close all;

addpath('../VEM2D/')

%%  PROBELM
a = 2;
f = @(X) (4*pi^2-a^2)*exp(-a*X(:,1)).*cos(2*pi*X(:,2));
gD = @(X) exp(-a*X(:,1)).*cos(2*pi*X(:,2));
gN = @(X) a*gD(X);

%%  GRID DIMENSIONS

nVec = [8, 16, 32, 64];
nIt = numel(nVec);
errVec = zeros(nIt, 3);

for i = 1:nIt

    %% GENERATE GRID
    
    n = nVec(i);
    G = unitSquare([n,n],[1,1]);
    G = sortEdges(G);
    G = computeVEM2DGeometry(G);

    %%  SET BC
    
    boundaryEdges = find(any(G.faces.neighbors == 0,2));
    tol = 1e-10;
    isNeu = abs(G.faces.centroids(boundaryEdges,1)) < tol;
    bc = VEM2D_addBC([], G, boundaryEdges(~isNeu), 'pressure', gD);
    bc = VEM2D_addBC(bc, G, boundaryEdges(isNeu), 'flux', gN);

    %%  SOLVE
    
    [sVEM1, G1] = VEM2D(G,f,1,bc, 'projectors', true);
    [sVEM2, G2] = VEM2D(G,f,2,bc, 'projectors', true);

    %%  CALUCLATE ERROR
    
    h = mean(G.cells.diameters);
    area = sqrt(sum(G.cells.volumes.^2));
    nK = G.cells.num;

    l2Err1 = l2Error(G1, sVEM1, gD, 1);
    l2Err2 = l2Error(G2, sVEM2, gD, 2);
    errVec(i,:) = [h, sqrt(sum(l2Err1)), sqrt(sum(l2Err2))];

    %%  PLOT GRID AND SOLUTIONS

    %   Grid
    
    gridFig = figure;
    set(gridFig, 'visible','off')
    plotGrid(G, 'facealpha', .2);
    set(gridFig,'DefaultTextInterpreter', 'LaTex');
    set(gca, 'XTick', [0,1]);
    set(gca, 'YTick', [0,1]);
    xlabel('$x$'); ylabel('$y$');
    axis equal off;
    
    cut = 4;
    ps = get(gridFig, 'Position');
    ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
    paperWidth = 10;
    paperHeight = paperWidth*ratio - cut;
    set(gridFig, 'paperunits', 'centimeters');
    set(gridFig, 'papersize', [paperWidth paperHeight]);
    set(gridFig, 'PaperPosition', [0    0   paperWidth paperHeight]);
    fileName = strcat('../../tex/thesis/fig/Grid_', num2str(i));
    print(gridFig, '-dpdf', fileName);

    %   1st order solution

    sol1Fig = figure;
    set(sol1Fig, 'visible','off')
    plotVEM2D(G, sVEM1, 1);
    set(sol1Fig,'DefaultTextInterpreter', 'LaTex');
    set(gca, 'XTick', [0,.5,1]);
    set(gca, 'YTick', [0,.5,1]);
    xlabel('$x$'); ylabel('$y$'); zlabel('$u_h$');
    zlim([-1.5,1.1]);
    
    cut = 4;
    ps = get(sol1Fig, 'Position');
    ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
    paperWidth = 10;
    paperHeight = paperWidth*ratio - cut;
    set(sol1Fig, 'paperunits', 'centimeters');
    set(sol1Fig, 'papersize', [paperWidth paperHeight]);
    set(sol1Fig, 'PaperPosition', [0    0   paperWidth paperHeight]);
    fileName = strcat('../../tex/thesis/fig/Sol1_', num2str(i));
    print(sol1Fig, '-dpdf', fileName);

    %   2nd order solution

    sol2Fig = figure;
    set(sol2Fig, 'visible','off')
    plotVEM2D(G, sVEM2, 2);
    set(sol2Fig,'DefaultTextInterpreter', 'LaTex');
    set(gca, 'XTick', [0,.5,1]);
    set(gca, 'YTick', [0,.5,1]);
    xlabel('$x$'); ylabel('$y$'); zlabel('$u_h$');
    zlim([-1.5,1.1]);
    
    cut = 4;
    ps = get(sol2Fig, 'Position');
    ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
    paperWidth = 10;
    paperHeight = paperWidth*ratio - cut;
    set(sol2Fig, 'paperunits', 'centimeters');
    set(sol2Fig, 'papersize', [paperWidth paperHeight]);
    set(sol2Fig, 'PaperPosition', [0    0   paperWidth paperHeight]);
    fileName = strcat('../../tex/thesis/fig/Sol2_', num2str(i));
    print(sol2Fig, '-dpdf', fileName);
    
    

end

%%  PLOT CONVERGENCE RATES

%   1st order convergence

conv1Fig = figure;
set(conv1Fig,'DefaultTextInterpreter', 'LaTex');
loglog(errVec(:,1), errVec(:,2), '-s');
hold on
loglog(errVec(:,1),(errVec(:,1)./errVec(end,1)).^2*errVec(end,2)*.8)
p1 = polyfit(log(errVec(:,1)), log(errVec(:,2)),1);
lStr = strcat('Slope = ', num2str(p1(1), '%.3f'));
h = legend(lStr, 'Slope =2.0');
set(h,'Interpreter','latex');
xlabel('$\log(h)$'); ylabel('$\log\left(\left\|u-\Pi^\nabla u_h\right\|_{0,\Omega}\right)$');
axis equal;

cut = 4;
ps = get(conv1Fig, 'Position');
ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
paperWidth = 10;
paperHeight = paperWidth*ratio - cut;
set(conv1Fig, 'paperunits', 'centimeters');
set(conv1Fig, 'papersize', [paperWidth paperHeight]);
set(conv1Fig, 'PaperPosition', [0    0   paperWidth paperHeight]);
fileName = strcat('../../tex/thesis/fig/Conv1');
print(conv1Fig, '-dpdf', fileName);

%   2nd order convergence

conv2Fig = figure;
set(conv2Fig,'DefaultTextInterpreter', 'LaTex');
loglog(errVec(:,1), errVec(:,3), '-s');
hold on
loglog(errVec(:,1),(errVec(:,1)./errVec(end,1)).^3*errVec(end,3)*.6)
p2 = polyfit(log(errVec(:,1)), log(errVec(:,3)),1);
lStr = strcat('Slope =', num2str(p2(1), '%.3f'));
h = legend(lStr, 'Slope =3.0');
set(h,'Interpreter','latex');
xlabel('$\log(h)$'); ylabel('$\log\left(\left\|u-\Pi^\nabla u_h\right\|_{0,\Omega}\right)$');
axis equal

cut = 4;
ps = get(conv2Fig, 'Position');
ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
paperWidth = 10;
paperHeight = paperWidth*ratio - cut;
set(conv2Fig, 'paperunits', 'centimeters');
set(conv2Fig, 'papersize', [paperWidth paperHeight]);
set(conv2Fig, 'PaperPosition', [0    0   paperWidth paperHeight]);
fileName = strcat('../../tex/thesis/fig/Conv2');
print(conv2Fig, '-dpdf', fileName);
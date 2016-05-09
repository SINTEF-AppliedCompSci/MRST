clc; clear all; close all;

addpath('../../../pebiGridding/voronoi3D/')
addpath('../VEM3D/')
addpath('../')

%%  PROBELM

f  = @(X) 4*pi^2*X(:,1).*sin(2*pi*X(:,2).*X(:,3)).*(X(:,2).^2 + X(:,3).^2);
gD = @(X) X(:,1).*sin(2*pi*X(:,2).*X(:,3));
gN = @(X) 2*pi*X(:,1).*X(:,2).*cos(2*pi*X(:,2).*X(:,3));

%%  GRID DIMENSIONS

nVec = [400, 800, 1600, 3200];
nIt = numel(nVec);
errVec = zeros(nIt, 3);

azel = [150,30];

for i = 1:nIt

    %% GENERATE GRID
    
    fprintf('Generating grid ...\n')
    tic;
    n = nVec(i);
    G = voronoiCube(n,[1,1,1]);
    fprintf('Done in %f seconds.\n\n', toc);
    
    G = computeVEM3DGeometry(G);
    
    %%  SET BC
    
    boundaryEdges = find(any(G.faces.neighbors == 0,2));
    tol = 1e-10;
    isNeu = abs(G.faces.centroids(boundaryEdges,3)-1) < tol;
    bc = VEM3D_addBC([], boundaryEdges(~isNeu), 'pressure', gD);
    bc = VEM3D_addBC(bc, boundaryEdges(isNeu) , 'flux'    , gN);

    %%  SOLVE
    


    %%  CALUCLATE ERROR
    
    h = max(G.cells.diameters);
    area = sqrt(sum(G.cells.volumes.^2));
    nK = G.cells.num;


    
    %%  PLOT GRID AND SOLUTIONS

    Kc = G.cells.centroids;
    cells = 1:G.cells.num;
    r = .8; c = [1,1,0];
    cells = cells(sum(bsxfun(@minus, Kc, c).^2,2) > r^2);
    faceNum = mcolon(G.cells.facePos(cells),G.cells.facePos(cells+1)-1);
    faces = G.cells.faces(faceNum);
    
    remCells = 1:G.cells.num;
    remCells = remCells(~ismember(remCells, cells));
    
    GP = removeCells(G,remCells);
    outerFaces = any(GP.faces.neighbors == 0,2);
    
    outerFaces = find(ismember(G.faces.centroids,GP.faces.centroids(outerFaces,:),'rows'));
    
    clear GP Kc;
    
    %   Grid
    
    gridFig = figure;
    set(gridFig, 'visible','off')
    plotGrid(G, cells, 'facecolor', [238,232,170]/255);
    set(gridFig,'DefaultTextInterpreter', 'LaTex');
    set(gca, 'XTick', [0,1]);
    set(gca, 'YTick', [0,1]);
    xlabel('$x$'); ylabel('$y$');
    view(azel)
    axis equal off;
    
    cut = 4;
    ps = get(gridFig, 'Position');
    ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
    paperWidth = 10;
    paperHeight = paperWidth*ratio - cut;
    set(gridFig, 'paperunits', 'centimeters');
    set(gridFig, 'papersize', [paperWidth paperHeight]);
    set(gridFig, 'PaperPosition', [0    0   paperWidth paperHeight]);
    fileName = strcat('../../tex/thesis/fig/Grid3D_', num2str(i));
    print(gridFig, '-dpdf', fileName);
    
    %   1st order solution
    
    [sVEM1, G] = VEM3D(G,f,bc,1,'cellProjectors', true);
    l2Err1 = l2Error3D(G, sVEM1, gD, 1);
    
    
    if i == 1 || i == nIt
        
        sol1Fig = figure;
        set(sol1Fig, 'visible','off')
        sVEM1 = calculateFaceAverages(G,sVEM1);
        plotFaces(G,outerFaces,sVEM1.faceMoments(outerFaces));
        set(sol1Fig, 'DefaultTextInterpreter', 'LaTex');
        colorbar;
        view(azel);
        axis equal;
        set(gca, 'XTick', [0,.5,1]);
        set(gca, 'YTick', [0,.5,1]);
        set(gca, 'ZTick', [0,.5,1]);
        xlabel('$x$'); ylabel('$y$'); zlabel('$z$');
        colorbar

        cut = 4;
        ps = get(sol1Fig, 'Position');
        ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
        paperWidth = 10;
        paperHeight = paperWidth*ratio - cut;
        set(sol1Fig, 'paperunits', 'centimeters');
        set(sol1Fig, 'papersize', [paperWidth paperHeight]);
        set(sol1Fig, 'PaperPosition', [0    0   paperWidth paperHeight]);
        fileName = strcat('../../tex/thesis/fig/Sol3D1_', num2str(i));
        print(sol1Fig, '-dpdf', fileName);

    end
    
    clear sVEM1;

    %   2nd order solution
     
    [sVEM2, G] = VEM3D(G,f,bc,2,'cellProjectors', true); 
    l2Err2 = l2Error3D(G, sVEM2, gD, 2);
    
    if i == 1 || i == nIt

        sol2Fig = figure;
        set(sol2Fig, 'visible','off')
        plotFaces(G,outerFaces,sVEM2.faceMoments(outerFaces));
        set(sol2Fig,'DefaultTextInterpreter', 'LaTex');
        colorbar;
        set(gca, 'XTick', [0,.5,1]);
        set(gca, 'YTick', [0,.5,1]);
        set(gca, 'ZTick', [0,.5,1]);
        xlabel('$x$'); ylabel('$y$'); zlabel('$z$');
        view(azel)
        axis equal
        colorbar;

        cut = 4;
        ps = get(sol2Fig, 'Position');
        ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
        paperWidth = 10;
        paperHeight = paperWidth*ratio - cut;
        set(sol2Fig, 'paperunits', 'centimeters');
        set(sol2Fig, 'papersize', [paperWidth paperHeight]);
        set(sol2Fig, 'PaperPosition', [0    0   paperWidth paperHeight]);
        fileName = strcat('../../tex/thesis/fig/Sol3D2_', num2str(i));
        print(sol2Fig, '-dpdf', fileName);

    end
    
    clear sVEM2;
    
    errVec(i,:) = [h, sqrt(sum(l2Err1)), sqrt(sum(l2Err2))];

    clear gridFig sol1Fig sol2Fig G l2Err1 l2Err2 isNeu cells boundaryEdges faceNum faces outerFaces remCells;
    close all;
    
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

%%

cut = 4;
ps = get(conv1Fig, 'Position');
ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
paperWidth = 10;
paperHeight = paperWidth*ratio - cut;
set(conv1Fig, 'paperunits', 'centimeters');
set(conv1Fig, 'papersize', [paperWidth paperHeight]);
set(conv1Fig, 'PaperPosition', [0    0   paperWidth paperHeight]);
fileName = strcat('../../tex/thesis/fig/Conv3D1');
print(conv1Fig, '-dpdf', fileName);

%%

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

%%
cut = 4;
ps = get(conv2Fig, 'Position');
ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
paperWidth = 10;
paperHeight = paperWidth*ratio - cut;
set(conv2Fig, 'paperunits', 'centimeters');
set(conv2Fig, 'papersize', [paperWidth paperHeight]);
set(conv2Fig, 'PaperPosition', [0    0   paperWidth paperHeight]);
fileName = strcat('../../tex/thesis/fig/Conv3D2');
print(conv2Fig, '-dpdf', fileName);
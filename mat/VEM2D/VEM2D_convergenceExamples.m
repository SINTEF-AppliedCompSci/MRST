%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------

clc; clear; close all;

ex = 3;

neuEx = false;
switch ex
    case 1
        girdType = 'cart';
        f = @(X) zeros(size(X,1),1);
        C = -[.5,.5];
        gD = @(X) -log(1./(sqrt(sum(bsxfun(@minus, X, C).^2,2))));
    case 2
        gridType = 'pebi';
        f = @(X) sin(X(:,1));
        gD = @(X) sin(X(:,1));
    case 3
        gridType = 'pebi';
        plotSol = true;
        neuEx = true;
        a = 2;
        f = @(X) (4*pi^2-a^2)*exp(-a*X(:,1)).*cos(2*pi*X(:,2));
        gD = @(X) exp(-a*X(:,1)).*cos(2*pi*X(:,2));
        gN = @(X) a*gD(X);
    case 4
        neuEx = true;
        f = @(X) -6*X(:,1);
        gD = @(X) X(:,1).^3;
        gN = @(X) -3*X(:,1).^2;
    case 5
        gridType = 'pebi';
        plotSol = true;
        neuEx = true;
        f = @(X) pi^2*X(:,1).*sin(pi*X(:,2));
        gD = @(X) X(:,1).*sin(pi*X(:,2));
        gN = @(X) -sin(pi*X(:,2));
    case 6
        gridType = 'pebi';
        plotSol = true;
        neuEx = true;
        a = 2;
        b = 2;
        c = 10;
        gD = @(X) cos(a*pi*X(:,1)).*sin(b*pi*X(:,2)).*(X(:,1)-.5).^2*c;
        f = @(X) cos(a*pi*X(:,1)).*sin(b*pi*X(:,2)).*...
                 (a^2*pi^2*c*(X(:,1)-.5).^2 + 4*a*pi*c*(X(:,1)-.5) -2*c + b^2*pi^2*c*(X(:,1)-.5).^2);
        gN = @(X) sin(b*pi*X(:,2)).*(a*pi*sin(a*pi*X(:,1))*c.*(X(:,1)-.5).^2 - 2*c*cos(a*pi*X(:,1)).*(X(:,1)-.5));
end


nVec = [10, 20, 40, 80];
nIt = numel(nVec);
errVec = zeros(nIt, 3);

for i = 1:nIt
    
    n = nVec(i);
    if strcmp(gridType,'cart')
        G = cartGrid([n,n],[1,1]);
    else
        G = unitSquare([n,n],[1,1]);
    end
    G = sortEdges(G);
    G = computeVEM2DGeometry(G);

    boundaryEdges = find(any(G.faces.neighbors == 0,2));
    if neuEx
        tol = 1e-10;
        isNeu = abs(G.faces.centroids(boundaryEdges,1)) < tol;
    else
        isNeu = false(numel(boundaryEdges),1);
        gN = gD;
    end
    bc = VEM2D_addBC([], G, boundaryEdges(~isNeu), 'pressure', gD);
    bc = VEM2D_addBC(bc, G, boundaryEdges(isNeu), 'flux', gN);
    
    [sVEM1, G1] = VEM2D(G,f,1,bc, 'projectors', true);
    [sVEM2, G2] = VEM2D(G,f,2,bc, 'projectors', true);
    
    h = mean(G.cells.diameters);
    area = sqrt(sum(G.cells.volumes.^2));
    nK = G.cells.num;
        
    l2Err1 = l2Error(G1, sVEM1, gD, 1);
    l2Err2 = l2Error(G2, sVEM2, gD, 2);
    errVec(i,:) = [h, sqrt(sum(l2Err1)), sqrt(sum(l2Err2))];
    
    if plotSol
        
        %   Grid
        
        gridFig = figure;
        plotGrid(G, 'facealpha', .2);
        set(gcf,'DefaultTextInterpreter', 'LaTex');
        set(gca, 'XTick', [0,1]);
        set(gca, 'YTick', [0,1]);
        xlabel('$x$');
        ylabel('$y$');
        axis([-.1, 1.1, -.1 1.1])
        axis equal;
        
        %   1st order solution
        
        sol1Fig = figure;
        plotVEM2D(G, sVEM1, 1);
        set(gcf,'DefaultTextInterpreter', 'LaTex');
        xlabel('$x$'); ylabel('$y$'); zlabel('$u_h$');
        
        %   2nd order solution
        
        sol2Fig = figure;
        plotVEM2D(G, sVEM2, 2);
        set(gcf,'DefaultTextInterpreter', 'LaTex');
        xlabel('$x$'); ylabel('$y$'); zlabel('$u_h$');
        
        pause;
        
    end
end

loglog(errVec(:,1), errVec(:,2), errVec(:,1), errVec(:,3));
legend('VEM 1st order', 'VEM 2nd order');
p1 = polyfit(log(errVec(:,1)), log(errVec(:,2)),1);
p2 = polyfit(log(errVec(:,1)), log(errVec(:,3)),1);

p1(1)
p2(1)
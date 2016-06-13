%   Poisson problem examples for VEM2D.
%
%   Example 2: 
%
%   Example 1: Point source at (0.5, 0.5)^T.
%
%-----------------------------------------------------------------ØSK-2016-

clc; clear; close all;

ex = 1;

switch ex
    
    case 1
                
        load('showInnerOuterCellsElephant.mat');
        G = sortEdges(G);
        G = computeVEM2DGeometry(G);
        
        f = @(X) X(:,1).*cos(5*pi*X(:,2));
        gD = @(X) 25*pi^2*X(:,1).*cos(X(:,2));
        
        bcEdges = find(any(Gs.faces.neighbors == 0, 2));
        bc = VEM2D_addBC([], G, bcEdges, 'pressure', gD);
        k = 2;
        
        sol = VEM2D(G,f,bc,k, 'src', src);
        plotCellData(G, sol.cellMoments);
        set(gcf, 'defaultTextInterpreter', 'LaTex');
        xlabel('$x$'); ylabel('$y$');
    
    case 3
        
        f = 0;
        C = [0.5, 0.5];
        gD = @(X) -log(sqrt(sum(bsxfun(@minus,X,C).^2,2)))/(2*pi)    ;
        gN = @(X) -(X(:,2)-C(2))./(2*pi*sum(bsxfun(@minus,X,C).^2,2));
        
        load singularityGrid.mat;      
        G = sortEdges(G);
        G = computeVEM2DGeometry(G);

        tol = 1e-6;
        xMax = 1; yMax = 1;

        bDir     = find(abs(G.faces.centroids(:,1))      < tol | ...
                        abs(G.faces.centroids(:,1)-xMax) < tol);
        bNeuS    = find(abs(G.faces.centroids(:,2))      < tol);
        bNeuN    = find(abs(G.faces.centroids(:,2)-yMax) < tol);

        bc = VEM2D_addBC([], G, bDir , 'pressure', gD         );
        bc = VEM2D_addBC(bc, G, bNeuS, 'flux'    , @(X) -gN(X));
        bc = VEM2D_addBC(bc, G, bNeuN, 'flux'    , gN         );
        
        Q = 1;
        srcCells = find(G.cells.tag);
        src = addSource([],srcCells(1),Q);
        
        k = 2;
        
        sol = VEM2D(G,f,bc,k, 'src', src);
        plotVEM2D(G,sol,k);
        axis([0,1,0,1]);
        set(gcf, 'defaultTextInterpreter', 'LaTex');
        xlabel('$x$'); ylabel('$y$'); zlabel('$u_h$');

        
    case 2
        xMax = 1; yMax = 1;
        f = @(X) sin(X(:,1));
        gD = @(X) sin(X(:,1));
    case 3
        neuEx = true;
        f = @(X) 24*exp(X(:,1)).*cos(5*X(:,2));
        gD = @(X) exp(X(:,1)).*cos(5*X(:,2));
        gN = @(X) -gD(X);
end


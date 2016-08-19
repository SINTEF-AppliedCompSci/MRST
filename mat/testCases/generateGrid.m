clc; clear; close all;

addpath('~/Documents/master/pebiGridding/voronoi2D/');
addpath('../VEM2D/');

%% Chose grid
% Chose the type of grid.
% gT = 1      Coarse cartesian
% gT = 2      composite pebi
% gT = 3      fully unstructured grid

xMax = 1; yMax = 1;

C = [xMax/2, yMax/2];
m = 10;
wellLine = {C};                % Set source center

gT = 3;

switch gT
    case 1

        nx = 61;
        G = cartGrid([nx,nx],[xMax,yMax]);
        G = computeGeometry(G);
        w1 = wellLine{1};
        D = pdist2(G.cells.centroids, w1);
        [~, I] = min(D, [], 1);
        G.cells.tag = false(G.cells.num,1);
        G.cells.tag(I') = true(size(I'));
        G = sortEdges(G);
        G = computeVEM2DGeometry(G);
        
        save('cartesianGrid.mat', 'G');

    case 2 

        gridSize = xMax/10;                   % Size of gridcells
        mlqtMax = 2;                            % Set number of reminement levels
        wellGridSize = 0.75/2^mlqtMax;          % Set relative well grid size
        mlqtSizes = 2.0*linspace(gridSize,gridSize*wellGridSize,mlqtMax+1)';
                                                % Size around wells to be refined 
        G = compositePebiGrid([gridSize, gridSize], [xMax, yMax], 'wellLines', wellLine, ...
                             'wellGridFactor', wellGridSize, ...
                             'mlqtMaxLevel', 2, 'mlqtLevelSteps', mlqtSizes);
        G = sortEdges(G);
        G = computeVEM2DGeometry(G);
        
        save('compPebiGrid.mat', 'G');
        
        plotGrid(G);

    case 3
        
        gridSize = xMax/10;
        wellGridSize = 0.7/2^2;
        epsilon = gridSize*.7;
        G = pebiGrid(gridSize, [xMax, yMax], 'wellLines', wellLine,   ...
                    'wellGridFactor', wellGridSize, 'wellRefinement',true);
        G = sortEdges(G);
        m = 1;
        G = computeVEM2DGeometry(G);

        save('untructPebi.mat','G');
        plotGrid(G);
        
    otherwise

    error('unknown grid case')

end      
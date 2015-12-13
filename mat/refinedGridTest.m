%% Grid with single source centered at 0.5,0.5
close all; clear;
%addpath('../../coop/Numerical-testing-of-new-features-in-MRST/')
path = fullfile('../../coop/Numerical-testing-of-new-features-in-MRST/','distmesh/');
mrstPath('reregister','distmesh',path);
mrstModule add distmesh

xmax = 1;                              % Set grid dimentions
ymax = 1; 
gridSize = 1/20;                       % Set size of grid cells
% delta = 0.001;
% gD = @(X) 10*ones(size(X,1),1);
% gN = @(X) zeros(size(X,1),1);
% f = @(X) 1/sqrt(2*delta*pi).*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.5).^2)./(2*delta));
% 
% boundaryEdges = find((G.faces.neighbors(:,1) == 0) + (G.faces.neighbors(:,2) == 0));
% bNeu = zeros(numel(boundaryEdges),1);
% for e = 1:numel(boundaryEdges)
%     X = G.nodes.coords(G.faces.nodes(G.faces.nodePos(e):G.faces.nodePos(e +1)-1));
%     if X(:,1) < 0.001
%         bNeu(e) = 1;
%     end
% end
% bc = struct('bcFunc', {{gN,gD}}, 'bcFaces', {{boundaryEdges(bNeu == 1), boundaryEdges(bNeu == 0)}}, 'bcType', {{'neu', 'dir'}});
faultLine = {};
wellLine = {[0.5,0.5]};                % Set source center

mlqtMax = 2;                        % Set number of reminement levels
wellGridSize = 0.5/2^mlqtMax;       % Set relative well grid size
mlqtSizes = 2*linspace(gridSize,gridSize*wellGridSize,mlqtMax+1)';
                                    % Set distance around wells to be
                                    % refined.

wellEps = 1/4;                         % Size around wells to be refined
                                       % (For unstructured grid)


% Create semi-structured grid
Gp = compositeGridPEBI(gridSize, [xmax, ymax], 'wellLines', wellLine, ...
                      'wellGridFactor', wellGridSize, ...
                      'mlqtMaxLevel', 2, 'mlqtLevelSteps', mlqtSizes);
                  
% Create fully unstructured grid
Gdist = compositeGridPEBIdistmesh(1/19, [1, 1], 'wellLines', wellLine, ...
                                 'wellGridFactor', 0.02*19, 'wellRefDist',wellEps);

run('../../matlab/project-mechanics-fractures/mystartup.m')
G = Gdist;
G = mrstGridWithFullMappings(G);
G = computeGeometry(G);

plotGrid(Gdist)

delta = 0.001;
gD = @(X) 10*ones(size(X,1),1);
f = @(X) -1/sqrt(2*delta*pi).*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.5).^2)./(2*delta));

tol = 1e-6;
boundaryEdges = find((G.faces.neighbors(:,1) == 0) + (G.faces.neighbors(:,2) == 0));
bc = struct('bcFunc', {{gD}}, 'bcFaces', {{boundaryEdges}}, 'bcType', {{'dir'}});

U = VEM2D(G,f,bc);

plotVEM(G, U, '')
%% Single fault
% close all; clear
% xmax = 1;                              % Set grid dimentions
% ymax = 1; 
% gridSize = 1/20;                       % Set size of grid cells
% 
% x1 = linspace(0.2,0.8,10);
% y1 = [0.2,0.25,0.3,0.34,0.40,0.5,0.6,0.65,0.7,0.85];
% faultLine = {[x1',y1']};                % Set fault line
% wellLine = {};
% 
% 
% faultGridSize = 0.5;                   % Set relative well grid size
% 
% % Create semi-structured grid
% Gp = compositeGridPEBI(gridSize, [xmax, ymax], 'faultLines', faultLine,...
%                        'faultGridFactor', faultGridSize);
%                   
% % Create fully unstructured grid
% Gdist = compositeGridPEBIdistmesh(1/19, [1, 1], 'faultLines', faultLine,...
%                                  'faultGridFactor', faultGridSize);

                             

                             
                             
                             
% %% Ploting
% figure()
% hold on
% plotGrid(Gp, 'faceColor', 'none')
% axis equal tight off
% hold on
% 
% %Plot cells taged as wells.
% if (numel(faultLine)>0), plotFault(Gp); end
% if (numel(wellLine)>0), plotWells(Gp);end
% 
% for i = 1:numel(wellLine)
%   line = wellLine{i};
%   if size(line,1) == 1
%       plot(line(1,1), line(1,2),'.r', 'markersize', 8);
%   end
%   plot(line(:, 1), line(:, 2),'r');
% end
% 
% 
% figure()
% hold on
% plotGrid(Gdist, 'faceColor', 'none')
% axis equal tight off
% hold on
% if numel(faultLine)>0, plotFault(Gp); end
% if numel(wellLine)>0, plotWells(Gp);end
% for i = 1:numel(wellLine)
%   line = wellLine{i};
%   if size(line,1) == 1
%       plot(line(1,1), line(1,2),'.r', 'markersize', 8);
%   end
%   plot(line(:, 1), line(:, 2),'r');
% end
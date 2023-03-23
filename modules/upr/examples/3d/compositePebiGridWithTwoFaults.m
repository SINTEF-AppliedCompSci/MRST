%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016-2020 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%% Example of two faults intersecting
% In this example, we create a grid using the wrapping function 
% compositePebiGrid3D.

%% Create Faults
% We first need to create a triangulation of the fault surfaces
nx = 10;
x = linspace(0,1,nx)';
z = x;
[X,Z] = meshgrid(x,z);
fDt1.ConnectivityList = delaunay([X(:),Z(:)]);
fDt2.ConnectivityList = delaunay([X(:),Z(:)]);

fDt1.Points = [X(:), 0.5*X(:)+0.1, Z(:)]*0.99 + 0.005;
fDt2.Points = [X(:),-0.5*X(:)+0.6, Z(:)];
rho = @(p) 1/nx*ones(size(p,1),1);

%% Create Grid
% We can now create the grid.
G = compositePebiGrid3D([nx,nx,nx],[1,1,1],'faceConstraints',{fDt1,fDt2}, ...
                                           'FCRho',{rho});
G = computeGeometry(G);
%% Plot Grid
% Let us plot the grid. Some cells are removed to better visualizing what
% the gridding of the faults.
figure()
cr = G.cells.centroids(:,2)<-0.5*G.cells.centroids(:,1) + 0.6 ...
    &G.cells.centroids(:,2)> 0.5*G.cells.centroids(:,1) + 0.1 ...
    &G.cells.centroids(:,3)<0.8;
plotGrid(G, ~cr)
axis equal tight
view(-83,50)
camlight

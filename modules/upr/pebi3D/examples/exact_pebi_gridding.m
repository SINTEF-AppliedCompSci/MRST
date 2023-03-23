% This script shows how to create a pebi grids of different dimensions
% where the higher dimensional grids conform exactly to the lower
% dimensional grids. The 3D gridding is not as automatic as the 2D
% gridding, but this script should show all the neccesary steps in creating
% the grids.


%% Define the fautls
% We start by defining the faults. The fractures are defined by a set of
% vertices. First we define two square faults:
f1 = [1,3,2; 4,3,2; 4,3,4; 1,3, 4];
f2 = [2,2,3.3; 5,2,3.3; 5,4,3.3; 2,4, 3.3];

% Then we use the help function ellipticSurface3D to create an elliptic
% fault
f3 = ellipticSurface3D([3,2.5, 3.2], 2, 2, 0, pi/2, pi/2);

fracs = {f1, f2, f3};

%% Finding fault intersections

intersections = surfaceIntersections3D(fracs);

%% Generate 1D grids
% We start by generating the 1D grids. The only parameter we have to
% specify here is the grid size of the 1D cells. The 1D grids are defined
% by the intersection of two faults. To generate the grids we call the 
% function lineGrid3D. This returns a cell array of all 1D grids.
ds_1 = 0.2;
grids_1 = lineGrid3D(intersections, ds_1);

%% Generate 2D grids
% We now generate the 2D grids. Here we specify two parameters. The first
% is the grid size of the 2D cells, the second is how far from the 1D lines
% the fault sites should be placed. If you have problems at the 
% intersection of two faults, try to reduce this parameter
ds_2 = 0.25;
gamma_2 = ds_2/6;

[grids_2] = surfaceGrid3D(fracs, grids_1, intersections, ds_2, gamma_2);

%% Plot grids so far.
% Let us plot the grids we have generated so far:
figure(1); clf; hold on
for i = 1:numel(grids_2)
   g = grids_2{i};
   plotGrid(g,'faceColor','r');
end
axis equal tight
%% Generate 3D grids
% Ok, so this looks good so far. Let us now generate the 3D grids
ds_3 = .4;
gamma_3 = ds_3/12;
grids_3 = volumeGrid3D([6,6,6], fracs, grids_2, ds_3, gamma_3);
grids_3 = computeGeometry(grids_3);

%% Plot Grids
c_to_plt = grids_3.cells.centroids(:,1) > 3;
plotGrid(grids_3, c_to_plt)

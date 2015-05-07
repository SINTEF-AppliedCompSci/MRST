%% Grids from standard data sets in MATLAB

load trimesh2d;
G = triangleGrid([xfe yfe], trife); figure, plotGrid(G);

load tetmesh;
G = tetrahedralGrid(X,tet); figure; plotGrid(G); view(3);
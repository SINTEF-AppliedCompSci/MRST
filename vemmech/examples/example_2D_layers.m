%% Example using a complex grid to test mechanical solver
%  It is possible to vary the material parametes in layers
%  
E = [1; 10; 1]; % define youngs modulo of layres
nu = [0.3; 0.3; 0.3]; % define Poison ratio of layers
zlayers = [0.3; 0.3; 0.3]; % values of layers
zlayers = 1/sum(zlayers)*zlayers; 
L = [1 1]; % define domain size
dims = [10 10]; % define grid dimensions

% Using the utiliy function squareLayerTests several types of grids are
% possible se this function for definitions

% define grids
% grid_type = 'cartgrid'; % define type of trid
% grid_type = 'triangle';
grid_type = 'mixed4';dims= [4 4]; %this needs small grid size  

% also define testcase type i.e. boundary conditions and forces
% test_name = 'original_2D_test';
test_name = 'hanging_rod_2D'

disturb = 0.042;% use twister to distrube grid

% define test cases by a utility function
[G, bc, test_cases] = squareLayersTest('E', E, 'nu', nu, 'zlayers', zlayers, 'cartDims', ...
                                       dims, 'L', L, 'grid_type', grid_type, 'disturb', ...
                                       disturb, 'make_testcases', false, 'test_name', ...
                                       test_name);

                                   G = computeGeometry(G);


%  get the test case paramters
testcase = test_cases{1};
% boundary condition for mechanics
el_bc = testcase.el_bc;
% load
load = testcase.load;
% material parameters in Voigts notation
C = testcase.C;
% Solve the linear elastic system using virtual element methods
[uVEM, extra] = VEM_linElast(G, C, el_bc, load);

% plot results
plotopts = {'EdgeAlpha', 0.0, 'EdgeColor', 'none'};
plotopts1 = {'EdgeAlpha', 0.1};
figure(1), clf
subplot(3, 1, 1)
plotNodeData(G, uVEM(:, 1), plotopts{:}); colorbar();
title('x displacement')
subplot(3, 1, 2)
plotNodeData(G, uVEM(:, 2), plotopts{:}); colorbar();
title('y displacement')
subplot(3, 1, 3)
% calculate div operator
vdiv = VEM2D_div(G);
mdiv = vdiv*reshape(uVEM', [], 1)./G.cells.volumes;
title('div')
plotCellDataDeformed(G, mdiv, uVEM, plotopts1{:}); colorbar();



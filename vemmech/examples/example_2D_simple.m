%% Example using a complex grid on a square domain test mechanical solver
%  
%  
E = 1; % Youngs modulo
nu = 0.3; %Poison ratio
L = [1 1]; % grid domain
dims = [128 128]/8; % cartesian grid dimension
%grid_type = 'cartgrid';
grid_type = 'mixed4';
test_name = 'original_2D_test'; % list of test names to be found in file
                                % squareTest. Note that CC schemes do not work for all
                                % type of boundary conditions.
disturb = 0.03; % disturbe the grid slightly

% define square test case
[G, bc, test_cases] = squareTest('E', E, 'nu', nu, 'cartDims', dims, ...
                                 'L', L, 'grid_type', grid_type, 'disturb', disturb, ...
                                 'make_testcases', false, 'test_name', test_name);

G = computeGeometry(G);
Ev = repmat(E, G.cells.num, 1); 
nuv = repmat(nu, G.cells.num, 1);

testcase = test_cases{1};
el_bc = testcase.el_bc;
load = testcase.load;
C = testcase.C;
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
%plotGridDeformed(G, uVEM, plotopts1{:}), colorbar(), 
vdiv = VEM2D_div(G);
mdiv = vdiv*reshape(uVEM', [], 1)./G.cells.volumes;
title('div')
plotCellDataDeformed(G, mdiv, uVEM, plotopts1{:}); colorbar();


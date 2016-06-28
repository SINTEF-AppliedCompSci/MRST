E = [1; 10; 1];
nu = 1e0*[0.3; 0.3; 0.3];
zlayers = [0.3; 0.3; 0.3];
zlayers = 1/sum(zlayers)*zlayers;
L = [1 1];
dims = [128 128];
grid_type = 'cartgrid';
% grid_type = 'triangle';
test_name = 'original_2D_test';
disturb = 0;

[G, bc, test_cases] = squareLayersTest('E', E, 'nu', nu, 'zlayers', zlayers, 'cartDims', ...
                                       dims, 'L', L, 'grid_type', grid_type, 'disturb', ...
                                       disturb, 'make_testcases', false, 'test_name', ...
                                       test_name);
G = computeGeometryCalc(G);
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
vdiv = VEM2D_div(G);
mdiv = vdiv*reshape(uVEM', [], 1)./G.cells.volumes;
title('div')
plotCellDataDeformed(G, mdiv, uVEM, plotopts1{:}); colorbar();



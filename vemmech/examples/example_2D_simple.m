E = 1;
nu = 0.3;
L = [1 1];
dims = [128 128]/8;
grid_type = 'cartgrid';
grid_type = 'mixed4';
test_name = 'original_2D_test'; % list of test names to be found in file
                                % squareTest. Note that CC schemes do not work for all
                                % type of boundary conditions.
disturb = 0.03;

[G, bc, test_cases] = squareTest('E', E, 'nu', nu, 'cartDims', dims, ...
                                 'L', L, 'grid_type', grid_type, 'disturb', disturb, ...
                                 'make_testcases', false, 'test_name', test_name);

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
%plotGridDeformed(G, uVEM, plotopts1{:}), colorbar(), 
vdiv = VEM2D_div(G);
mdiv = vdiv*reshape(uVEM', [], 1)./G.cells.volumes;
title('div')
plotCellDataDeformed(G, mdiv, uVEM, plotopts1{:}); colorbar();

return

[uCC, out] = CC_linElast(G, C, el_bc, load);

% plot results
figure(2)
title('CC')
subplot(3, 1, 1)
plotCellData(G, uCC(:, 1)), colorbar(), 
subplot(3, 1, 2)
plotCellData(G, uCC(:, 2)), colorbar(), 

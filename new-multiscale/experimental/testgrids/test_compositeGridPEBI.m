%% Show 3d version
Gp = compositeGridPEBI([15, 15], [10, 10]);

figure;
plotGrid(Gp)

figure;
plotGrid(makeLayeredGrid(Gp, 10))
%% Local well refinement
close all

Wpt = [3, 3; 7, 16];

Gp = compositeGridPEBI([15, 30], [10, 20], 'Wells', Wpt, 'radnum', 5, ...
                                           'radius', 2.0, 'growthfactor', 1.3);
figure;
plotGrid(Gp)
axis equal tight off

%% Check geometry
figure;
Gp = computeGeometry(Gp);
plotCellData(Gp, Gp.cells.volumes);
colormap jet
axis equal tight off

%% Fault honoring pebi / triangle grids
close all
l = [0.50, 0.50; 0.8, 1.1];

Gp = compositeGridPEBI([15, 30], [1, 2], 'lines', {l}, 'makePEBI', true);
Gt = compositeGridPEBI([15, 30], [1, 2], 'lines', {l}, 'makePEBI', false);

figure;
plotGrid(Gp)
axis equal tight off
hold on
plot(l(:, 1), l(:, 2));


figure;
plotGrid(Gt)
axis equal tight off
hold on
plot(l(:, 1), l(:, 2));
%% Extra points locally
makepts = @(N, pos, dims) bsxfun(@plus, bsxfun(@times, -0.5 + rand(N, 2), dims), pos);
pts = [makepts(50, [5, 5], [2, 2]); makepts(25, [5, 5], [5, 5])];

Gp = compositeGridPEBI([15, 30], [10, 20], 'extraPts', pts);

figure;
plotGrid(Gp)
axis equal tight off
%% Multiple faults

% close all
l = [0.50, 0.50; 0.8, 1.1];
l2 = [0.1, 1.5; 0.4, 1.1];
faults = {l, l2};
Gp = compositeGridPEBI([15, 30], [1, 2], 'lines', faults, 'makePEBI', true);

figure;
plotGrid(Gp)
axis equal tight

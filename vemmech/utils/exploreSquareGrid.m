%% Explore the different types of grids that can be set up using the function squareGrid
%

figure()
plotGrid(squareGrid([4, 4, 3], [3, 3, 2], 'grid_type', 'cartgrid', 'disturb', ...
                    0.1));
title('cartgrid');
view([30, 30])

figure()
plotGrid(squareGrid([4, 4], [3, 3], 'grid_type', 'triangle', 'disturb', 0.1));
title('triangle');

figure()
plotGrid(squareGrid([4, 4], [3, 3], 'grid_type', 'pebi', 'disturb', 0.1));
title('pebi');

grid_types = {'boxes1', 'boxes2', 'boxes3', 'boxes4'};
for i = 1 : 4
    figure()
    plotGrid(squareGrid([4, 4], [3, 3], 'grid_type', grid_types{i}, 'disturb', ...
                        0.1));
    title(grid_types{i});
end

grid_types = {'mixed1', 'mixed2', 'mixed3', 'mixed4'};
for i = 1 : 4
    figure()
    plotGrid(squareGrid([4, 4], [3, 3], 'grid_type', grid_types{i}, 'disturb', ...
                        0.1));
    title(grid_types{i});
end

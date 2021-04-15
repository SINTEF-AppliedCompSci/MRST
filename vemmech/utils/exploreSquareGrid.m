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

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>

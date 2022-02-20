%% 3d plotting
pdims = [1 1 2*pi];
G = cartGrid([10, 10, 15], pdims);
G = computeGeometry(G);


perfCells = zeros(G.cartDims(3), 1);
% Make a spiraling well
[ii, jj, kk] = gridLogicalIndices(G);
for i = 1:G.cartDims(3)
    cells = find(kk == i);
    
    pts = G.cells.centroids(cells, 1:2);
    z = G.cells.centroids(cells, 3);
    
    pt = 0.25*[sin(z), cos(z)] + repmat([0.5, 0.5], size(z, 1), 1);
    
    dist = bsxfun(@minus, pt, pts);
    dist = sqrt(sum(dist.^2, 2));
    
    [~, ind] = min(dist);
    
    perfCells(i) = cells(ind);
end

rock = struct('perm', ones(G.cells.num, 1));
W = addWell([], G, rock, perfCells, 'radius', 0.001);

% Plot the well
figure;
for i = 1:4
    subplot(2, 2, i)
    plotGrid(G, 'facec', 'none', 'edgea', .1)
    view(-90, 0)
    axis tight off
    
    switch i
        case 1
            plotGrid(G, perfCells);
            title('Just cells')
        case 2
            plotWell(G, W)
            title('Old plotting')
        case 3
            plotWellData(G, W)
            title('New plotting')
        case 4
            plotWellData(G, W, {cos(G.cells.centroids(perfCells, 3)*4 + pi/4 )}, 'radialData', true, 'radius', 2)
            title('Plotting colorized by data + radial modifier')
    end
end



%% 2d plotting
G = cartGrid([10, 10]);
G = computeGeometry(G);

rock.perm = ones(G.cells.num, 1);

W = [];
W = addWell(W, G, rock, 1);
W = addWell(W, G, rock, 10);
W = addWell(W, G, rock, 91);
W = addWell(W, G, rock, 100);


figure;
subplot(1, 2, 1)
plotWellData(G, W, 'TextColor', 'w')
plotGrid(G, 'FaceAlpha', 0)
axis tight off

subplot(1, 2, 2)
plotWellData(G, W, {1, 5, 3, 2}, 'TextColor', 'w')
plotGrid(G, 'FaceAlpha', 0)
axis tight off

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

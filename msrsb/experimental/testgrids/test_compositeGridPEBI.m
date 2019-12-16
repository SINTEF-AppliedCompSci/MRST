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

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.
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

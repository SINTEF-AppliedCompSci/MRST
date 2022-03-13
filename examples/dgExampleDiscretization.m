%% Components of a Discontinuous Galerkin Discretization in MRST
% This educational example goes through the main components of a
% discontinuous Galerkin (dG) discretization in MRST
mrstModule add dg
saveeps = @(a,b) disp(b);  % Dummy function
savepng = @(a,b) disp(b);  % Dummy function

%% Make grid
% We construct a small PEBI grid
mrstModule add upr            % For generating PEBI grids
G = pebiGrid2D(1/8, [1,1]);   % PEBI grid
G = computeGeometry(G);       % Compute geometry
G = computeCellDimensions(G); % Compute cell dimensions

gray    = [1,1,1]*0.8; % For plotting
addText = @(x, string) text(x(:,1), x(:,2), string, 'fontSize', 12, 'HorizontalAlignment', 'left');

%% Plot the grid
% We plot the grid and look at some geometric properties of a single cell
c = 22; % Cell number 17
% Plot cell 17 and its topological neighbors
figure();
plotGrid(G, 'facecolor', 'none');
cn = G.faces.neighbors(any(G.faces.neighbors == c,2),:);
cn = cn(cn>0);
plotGrid(G, cn(:), 'facecolor', 'none', 'linewidth', 2); axis equal off
plotGrid(G, c, 'facecolor', gray, 'linewidth', 2);
drawnow(), pause(0.2), saveeps('discretization', 'grid');
% By calling omputeCellDimensions(G), we computed the dimensions of a the
% smallest bounding box aligned with the coordinate axes so that the
% centroid of the box coincides with the cell centroid. We will use this to
% define our basis functions
xc = G.cells.centroids(c,:); % Centroid of cell 10
dx = G.cells.dx(c,:);        % Bounding box dimensions of cell 10
figure()
hold on
plotGrid(G, c, 'facecolor', gray);
rectangle('Position', [xc - dx/2, dx]);
plot(xc(1), xc(2), '.k', 'markerSize', 30);
plot([xc(1); xc(1)], [xc(2) - dx(2)/2; xc(2) + dx(2)/2], '--k');
plot([xc(1) - dx(1)/2; xc(1) + dx(1)/2], [xc(2); xc(2)], '--k');
axis equal off
drawnow(), pause(0.2), saveeps('discretization', 'cell');

%% Construct discretization
mrstModule add dg                        % Load module
G    = computeCellDimensions(G);         % Compute cell dimensions
disc = DGDiscretization(G, 'degree', 2); % Construct dG(2) discretization

%% Inspect basis
% The basis functions are conveniently constructed by a Polynomial class. A
% polynomial is then described as a set of exponents k and weights w, and
% the class has overloaded operators for elementary arithmetic operations
% (addition, subtraction, multiplication, etc.), as well as tensor
% products, derivatives and gradients. We have construct a Legendre-type
% basis of polynomials of order <= 2.
disp(disc.basis)        % Inspect basis
disp(disc.basis.psi{3}) % The third basis function is simply y

%% Plot polynomial basis
% The basis functions are evaluated in local coordinates in each cell:
% $\psi(\hat{x}) = \psi(\frac{x - x_c}{\Delta_x/2})
h = plotDGBasis(G, disc.basis, c, 'edgealpha', .05);
for i = 1:numel(h)
    set(0, 'CurrentFigure', h(i));
    ax = gca;
    [ax.XTickLabel, ax.YTickLabel] = deal({});
    ax.FontSize = 20; caxis([-1 1.2])
    savepng('discretization', ['basis-', num2str(i)]);
end

%% Construct triangle cubature
% DG involves excessive evaluation of integrals. This is handeled by
% cubatures, which are implemented in the Cubature class. The simplest
% possible approach for polygonal cells is to subdivide the cells into
% simplices, and use a cubature rule for each simplex. This is implemented
% in the TriangleCUbature/TetrahedronCubature classes, which is the default
% cubature of the DGDiscretization
triCubature = disc.cellCubature;
figure(), plotCubature(G, triCubature, c, 'faceColor', gray); axis equal tight
savepng('discretization', 'cubature-tri');

%% Test the triangle cubature
[W, x, ~, cells] = triCubature.getCubature(c, 'cell'); % Get points/weights
x = triCubature.transformCoords(x, cells); % Transform to reference coords
% The integral is easily evalauted using the matrix W, which has the
% cubature weights w placed so that int(psi)dx = W*psi(x). The weights sum
% to one, so that by construction, the integral of the first basis function
% (which is constant) should equal one, whereas the linear basis functions
% (1 and 2) should be zero
cellfun(@(psi) W*psi(x), disc.basis.psi)

%% Construct moment fitting cubature
% The cubature above uses much more points than strictly needed. We can
% construct a more efficient cubature by moment fitting. In this approach,
% we use a known cubature to compute the integrals for all functions in the
% basis for the polyomials we want our rule to be exact for. Then, given an
% initial set of cubature points, we compute corresponding weights by
% linear least squares, successively removing points until it is no longer
% possible to take away more points. In this case, the result is a cubature
% with six points, which equals the number of basis functions
% Moment-fitting without redcution
mrstModule add vem vemmech
degree = 2;
mfCubature_nr = MomentFitting2DCubature(G, 2*degree, 'chunkSize', 1, 'reduce', false);
figure(), plotCubature(G, mfCubature_nr, c, 'faceColor', gray); axis equal tight
savepng('discretization', 'cubature-mf-nr');
% Moment-fitting with redcution
mfCubature = MomentFitting2DCubature(G, 2*degree, 'chunkSize', 1);
figure(), plotCubature(G, mfCubature, c, 'faceColor', gray); axis equal tight
savepng('discretization', 'cubature-mf');

%% Test the moment fitting cubature
[W, x, w, cells] = mfCubature.getCubature(c, 'cell');
x = mfCubature.transformCoords(x, cells);
cellfun(@(psi) W*psi(x), disc.basis.psi)

%% Construct line cubature
% We also need a cubature for integrals over cell faces. The cubature for a
% an interface is used for both cells sharing the interface, but the
% functions are evaluated from different sides depending on which cell we
% are considering
lineCubature = disc.faceCubature;
faces     = G.cells.faces(G.cells.facePos(c):G.cells.facePos(c+1)-1);
figure(), plotCubature(G, lineCubature, faces(3), 'smax', 10, ...
            'faceColor', gray, 'plotBoundingBox', false); axis equal tight
savepng('discretization', 'cubature-line');

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
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

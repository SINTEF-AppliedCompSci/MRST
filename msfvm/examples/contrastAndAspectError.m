%% Example demonstrating the effects of aspect ratios and contrast on msfvm
% The Multiscale Finite Volume method is based on approximating boundary
% conditions using lower dimensional solutions. This approach is local in
% nature and gives good results for relatively smooth heterogeneties, but
% is sensitive to high contrast, high aspect ratios.
%
% In this example we will show how a single dataset, realized with
% different levels of contrast and aspect ratios can demonstrate the limits
% of the method.
%
% Note that this example is intentionally difficult, placing the dual grid
% boundary over a high contrast region. This makes the boundary conditions
% resulting from the lower dimensional numerical solution a poor
% approximation of the global field, giving both a degenerate coarse scale
% operator and a poor prolongation operator if the contrast is large.
%

mrstModule add mrst-gui coarsegrid msfvm incomp

G = cartGrid([30 30]);
G = computeGeometry(G);

Cdims = [3,3];
p = partitionUI(G, Cdims);

CG = generateCoarseGrid(G, p);
CG = coarsenGeometry(CG);
DG = partitionUIdual(CG, Cdims);

% Read in an image for permeability
im = double(imread('perm_mono.gif'))';

clf;
plotCellData(G, im(:), 'EdgeColor', 'none')
view(0, 90)
axis tight
title('Permeability dataset with coarse grids')

% Plot dual coarse grid
plotGrid(G, DG.ee, 'FaceColor', 'none', 'EdgeColor', 'w')
% Plot primal coarse grid
outlineCoarseGrid(G, p, [1 1 1].*.6)

%% Define a wide range of numerical experiments
% We define a matrix of different contrasts and aspect ratios to get data
% to show how the solution quality varies.
%
% At each combination, we store the error 
%
% $$ e_{ms} =  \frac{\| p_{ms} - p_{ref} \|_2}{\|p_{ref}\|_2} $$
%
% for further use with plotting.

cexp = 0:.5:8;
contrasts = 10.^-(cexp);
aspects = [1 5 10 25 100];
nc = numel(contrasts);
na = numel(aspects);


solutions = cell(nc, na);
errors = nan(nc, na);

normerr = @(ms, ref, k) norm(ref - ms, k)/norm(ref, k);

for i = 1:nc
    for j = 1:na
        aspect = aspects(j);
        % Set up grid with aspect ratio
        G = cartGrid([30 30], [30*aspect, 30]*meter);
        G = computeGeometry(G);
        
        % Scale permeability mapping with contrast ratio
        hi = 1*darcy;
        lo = contrasts(i)*darcy;

        im = double(imread('perm_mono.gif'))';
        im = im./max(im(:));

        im = hi*im + lo;
        rock.perm = double(im(:));

        T = computeTrans(G, rock);
        
        % Set up linear flow boundary conditions
        bc = [];
        d  = abs(G.faces.centroids(:,1) - max(G.faces.centroids(:,1)));
        es = find (d < 1e6*eps);

        d  = abs(G.faces.centroids(:,1) - min(G.faces.centroids(:,1)));
        ws = find (d < 1e6*eps);

        bc = addBC(bc, ws, 'pressure', 1);
        bc = addBC(bc, es, 'pressure', 0);
        
        % Trivial fluid
        fluid        = initSingleFluid('mu' ,    1*centi*poise     , ...
                                       'rho', 1014*kilogram/meter^3);

        sol = initResSol(G, 0);
        % Solve multiscale
        ms = solveMSFV_TPFA_Incomp(sol, G, CG, T, fluid, 'Dual', DG, 'Reconstruct', true, 'bc', bc);
        % Solve TPFA reference
        ref = incompTPFA(sol, G, T, fluid, 'bc', bc);
        % Store the error and solutions
        errors(i, j) = normerr(ref.pressure, ms.pressure, 2);
        solutions{i, j} = struct('ms', ms.pressure, 'ref', ref.pressure);
    end
end
%% Plot the numerical results
% We plot the error as a function of the increasing contrast for the
% different aspect ratios. As we can see, increasing either will make the
% problem more challenging.
%

clf

axis tight
plot(cexp, errors, '--', 'LineWidth', 2)
ylim([0 3])
legend(arrayfun(@(x) sprintf('Aspect ratio %d', x), aspects, 'unif', false), 'location', 'northwest')
xlabel('Contrast')
ylabel('|e|_2')


tickx = cexp(1:2:end);

set(gca, 'XTick', tickx, 'XTickLabel', arrayfun(@(x) sprintf('1e-%1.1g', x), tickx, 'unif', false));
axis tight
%% Plot the pressures
% We plot a subset of the solutions for both multiscale and reference, with
% the error as as seperate plot. Note that this forms a matrix of plots,
% where going from left to right along the rows or down the columns
% increases the error.
close all

h1 = figure('units','normalized','outerposition',[.5 .5 .5 .5]);
h2 = figure('units','normalized','outerposition',[.5 .5 .5 .5]);
h3 = figure('units','normalized','outerposition',[.5 .5 .5 .5]);

% csubs = 1:8:nc;
csubs = [3 5 7];
asubs = 1:2:na;

for i = 1:numel(csubs)
    for j = 1:numel(asubs)
        ci = csubs(i);
        aj = asubs(j);
        
        ms = solutions{ci, aj}.ms;
        ref = solutions{ci, aj}.ref;
        
        % Plot multiscale pressure
        set(0, 'CurrentFigure', h1);
        subplot(numel(csubs), numel(asubs), i + (j-1)*numel(csubs))
        plotCellData(G, ms, 'EdgeColor', 'none')
        axis tight off
        title(['MS pressure: Contrast ', num2str(contrasts(ci)), ' aspect ', num2str(aspects(aj)), ' error: ', num2str(errors(ci, aj))])
        
        caxis([min(ref), max(ref)])
        
        set(0, 'CurrentFigure', h2);
        subplot(numel(csubs), numel(asubs), i + (j-1)*numel(csubs))
        plotCellData(G, ref, 'EdgeColor', 'none')
        axis tight off
        title(['Reference pressure: Contrast ', num2str(contrasts(ci)), ' aspect ', num2str(aspects(aj))])

        set(0, 'CurrentFigure', h3);
        subplot(numel(csubs), numel(asubs), i + (j-1)*numel(csubs))
        plotCellData(G, abs(ms - ref)./ref, 'EdgeColor', 'none')
        axis tight off
        colorbar
        title(['Error in multiscale: Contrast ', num2str(contrasts(ci)), ' aspect ', num2str(aspects(aj)), ' error: ', num2str(errors(ci, aj))])
    end
end

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

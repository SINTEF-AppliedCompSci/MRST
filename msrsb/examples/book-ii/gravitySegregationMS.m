%% Simulate gravity segregation with the MsRSB method
% This example is a simple demonstration in how to use the multiscale
% solver to solve a gravity segregation problem using the MsRSB method.

mrstModule add coarsegrid incomp msrsb

%%  Set up the fine-scale model
rng(07521362);
gravity reset on
G    = computeGeometry(cartGrid([20 1 40],[500 20 60]*meter));
p    = gaussianField(G.cartDims, [0.2 0.4], [11 1 3], 2.5);
K    = p.^3.*(1e-5)^2./(0.81*72*(1-p).^2);
rock = makeRock(G, K(:), p(:));
hT   = computeTrans(G, rock);

figure
K = convertTo(rock.perm(:,1),milli*darcy);
plotCellData(G, K,'EdgeAlpha',.1); mrstColorbar(K,'south');
axis tight, view(0,0)

fluid = initSimpleFluid('mu', [1, .5]*centi*poise, ...
    'n', [2, 2], 'rho', [1000, 600]*kilo*gram/meter.^3);

state0 = initResSol(G, 0, [0, 1]);
ind    = G.cells.centroids(:,3)<30;
state0.s(ind,:) = 1 - state0.s(ind,:);


%% Set up multiscale problem
coarsen = [5 1 5];
coarsedims = ceil(G.cartDims./coarsen);
p     = partitionUI(G, ceil(G.cartDims./coarsen));
CG    = generateCoarseGrid(G, compressPartition(p));
CG    = coarsenGeometry(CG);
CG    = storeInteractionRegionCart(CG);
basis = getMultiscaleBasis(CG, getIncomp1PhMatrix(G, hT), 'type', 'msrsb');
smoother = getSmootherFunction('type', 'ilu0','iterations',1);
iargs = {'tolerance', 1e-3, 'iterations', 50,...
          'useGMRES', true, 'getSmoother', smoother};

%% Set up solvers
time = 10*year;
Nt   = 100;
dt   = time/Nt;

psolver  = @(state) incompTPFA(state, G, hT, fluid);
tsolver  = @(state) implicitTransport(state, G, dt, rock, fluid);
mssolve  = @(state) incompMultiscale(state, CG, hT, fluid, basis);
isolver  = @(state) incompMultiscale(state, CG, hT, fluid, basis, iargs{:});

[rstates,mstates,istates] = deal(cell(Nt+1,1));
rstates{1} = psolver(state0);
mstates{1} = mssolve(state0);
istates{1} = isolver(state0);

%% Plot the initial saturation distribution
figure('position', [50, 50, 1000, 500])
subplot(1,3,1)
hp(1) = plotCellData(G, rstates{1}.s(:,1),'EdgeAlpha',.1); 
view(0,0); axis tight, caxis([0 1]), colormap(flipud(winter))
title('Fine scale','FontSize',12,'FontWeight','normal');

subplot(1,3,2)
hp(2) = plotCellData(G, mstates{1}.s(:,1),'EdgeAlpha',.1); 
view(0,0); axis tight, caxis([0 1]), colormap(flipud(winter(10).^1.5))
plotFaces(CG,1:CG.faces.num,'EdgeColor','k','FaceColor','none');
title('MsRSB: no iterations','FontSize',12,'FontWeight','normal');

subplot(1,3,3)
hp(3) = plotCellData(G, istates{1}.s(:,1),'EdgeAlpha',.1); 
view(0,0); axis tight, caxis([0 1]), colormap(flipud(winter(10).^1.5))
plotFaces(CG,1:CG.faces.num,'EdgeColor','k','FaceColor','none');
title('MsRSB: iterative GMRES','FontSize',12,'FontWeight','normal');
[~,c] = boundaryFaces(G);

%% Run the simulation
% To accelerate rendering of new results, we do not replot the grid cells,
% but simply update the CData (saturation values) that are used to
% determine colors in each subplot.
for i = 1:Nt
    state        = psolver(rstates{i});
    rstates{i+1} = tsolver(state);
    hp(1).CData  = rstates{i+1}.s(c,1);
    
    state        = mssolve(mstates{i});
    mstates{i+1} = tsolver(state);
    hp(2).CData  = mstates{i+1}.s(c,1);
    
    state        = isolver(istates{i});
    istates{i+1} = tsolver(state);
    hp(3).CData  = istates{i+1}.s(c,1);
    drawnow
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

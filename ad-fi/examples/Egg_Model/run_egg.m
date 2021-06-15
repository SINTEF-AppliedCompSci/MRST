%% The Delft EGG model
% Researchers from TU Delft and Endhoven have developed their own model - the Egg model -
% to compare different two phase flow simulators. In this process, MRST has been validated
% against Eclipse, GPRS, and MoReS. The results obtained with the four simulators are
% almost identical. A detailed description of this project and its outcome can be found in
% the following <http://dx.doi.org/10.4121/uuid:916c86cd-3558-4672-829a-105c62985ab2
% research note>.



%% Read and process input files
% The data file for the egg model are available at
% <http://dx.doi.org/10.4121/uuid:916c86cd-3558-4672-829a-105c62985ab2
% 3TU.Datacentrum>. From there you can download the whole dataset but we will only need
% the following files which can be found under the directory
% |Egg_Model_Data_Files_v2/MRST/|:
% 
%  ACTIVE.INC,COMPDAT.INC, Egg_Model_ECL.DATA, mDARCY.INC, SCHEDULE_NEW.INC 
% 
% Put all these files under the same directory and let |fn| denote the path of
% |Egg_Model_ECL.DATA|. The next two lines correspond to the setup where the whole
% dataset has been downloaded under the same directory as |run_egg.m|.

require ad-fi deckformat ad-props

currentdir = fileparts(mfilename('fullpath'));
fn = fullfile(currentdir, 'data', 'Egg_Model_Data_Files_v2', 'MRST', 'Egg_Model_ECL.DATA');

deck = readEclipseDeck(fn);

% The deck is given in field units, MRST uses metric.
deck = convertDeckUnits(deck);

G = initEclipseGrid(deck);


%%
% The egg-shaped grid is constructed be removing cells from a rectangular cartesian
% grid. We use the MRST function extractSubgrid to select the active cells given in
% |deck.GRID.ACTNUM|. The active cells are loaded in the deck from |ACTIVE.INC|.

G = extractSubgrid(G, logical(deck.GRID.ACTNUM));
G = computeGeometry(G);

rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);

% Create a special ADI fluid which can produce differentiated fluid properties.
fluid = initDeckADIFluid(deck);

% The case includes gravity
gravity on


%% Plot wells and permeability
%

figure(1)
clf;
W = processWells(G, rock, deck.SCHEDULE.control(1));
plotCellData(G, convertTo(rock.perm(:,1), milli*darcy), 'FaceAlpha', .5, ...
            'EdgeAlpha', .3, 'EdgeColor', 'k');
plotWell(G, W, 'fontsize', 10, 'linewidth', 1);
title('Permeability (mD)')
axis tight;
view(35, 40);
colorbar('SouthOutside');

%% Set up the initial state
% We consider a reservoir with a uniform saturation distribution (s_w =
% 0.1, s_o = 0.9) and we compute an approximation of the initial pressure.

sw = 0.1; 
so = 0.9;
pr   = 400*barsa;
rz   = G.cells.centroids(1,3);
dz   = G.cells.centroids(:,3) - rz;
rhoO    = fluid.bO(400*barsa)*fluid.rhoOS;
rhoW    = fluid.bW(400*barsa)*fluid.rhoWS;
rhoMix  = sw*rhoW + so*rhoO;
p0   = pr + norm(gravity)*rhoMix*dz;

% The function initResSol initializes the solution data structure.
rSol  = initResSol(G, p0, [sw, so]);

system = initADISystem(deck, G, rock, fluid, 'cpr', true);
system.pscale = 1/(100*barsa);
system.nonlinear.cprBlockInvert = false;
system.nonlinear.cprRelTol      = 2e-2;
schedule = deck.SCHEDULE;


%% Run the simulation
% We use a fully implicit oil/water solver.

tt = tic;
[wellSols, states, iter] = runScheduleADI(rSol, G, rock, system, schedule);
toc(tt)


%% Plot oil and water production rates
%

% Put the well solution data into a format more suitable for plotting
[wrt, ort, grt, bhp] = wellSolToVector(wellSols);
T = convertTo(cumsum(deck.SCHEDULE.step.val), day);

figure(2)
clf
prodind = 9:12; % The four producers correspond to the four last entries in wrt, ort, grt.
for i = 1 : 4
   subplot(2,2,i);
   h = plot(T, ort(:, prodind(i))*day, T, wrt(:, prodind(i))*day);
   title(W(prodind(i)).name);
   legend({'Oil rate', 'Water rate'})
end

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

%% Add required modules
mrstModule add ad-core ad-blackoil ad-eor ad-props deckformat sequential

%% Setup example
fn    = 'POLYMER.DATA';
deck  = readEclipseDeck(fn);
deck  = convertDeckUnits(deck);
G     = initEclipseGrid(deck);
G     = computeGeometry(G);
rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);
fluid = initDeckADIFluid(deck);

% Create models
modelTFI = TwoPhaseOilWaterModel(G, rock, fluid, 'inputdata', deck);
modelPFI = OilWaterPolymerModel(G, rock, fluid, 'inputdata', deck);
modelTSQ = getSequentialModelFromFI(modelTFI);
modelPSQ = getSequentialModelFromFI(modelPFI);


%%

% Create schedule
schedule = convertDeckScheduleToMRST(modelPFI, deck);

% Reduce schedule
nsteps = 30;
scheduleOW = schedule;
scheduleOW.step.control  = 2.*ones(nsteps,1); % control 2 has no polymer
scheduleOW.step.val      = 10*day.*ones(nsteps,1);
scheduleOW.step.val(1)   = 3*day; % step
scheduleOW.step.val(2:3) = 5*day;

% Polymer schedule
scheduleP = scheduleOW;
scheduleP.step.control  = 1.*ones(nsteps,1); % control 1 has polymer

% Oil rel-perm from 2p OW system.
% Needed by equation implementation function 'eqsfiOWExplictWells'.
fluid.krO = fluid.krOW;

gravity reset on


%% Set up simulation parameters
% We want a layer of oil on top of the reservoir and water on the bottom.
% To do this, we alter the initial state based on the logical height of
% each cell. The resulting oil concentration is then plotted.
ijk = gridLogicalIndices(G);
state0 = initResSol(G, 300*barsa, [ .9, .1]);
state0.s(ijk{3} == 1, 2) = .9;
state0.s(ijk{3} == 2, 2) = .8;
state0.s(:,1) = 1 - state0.s(:,2); % Enforce s_w + s_o = 1;

% Add zero polymer concentration to the state.
state0.cp    = zeros(G.cells.num, 1);

%% Plot grid
figure;
plotCellData(G, rock.perm(:,1)./(milli*darcy), 'facealpha', 0.7);
view([-16,20]); axis tight; colorbar;
plotWell(G, scheduleOW.control(1).W);
xlabel('x-axis');ylabel('y-axis');zlabel('z-axis');


%% Run fully implicit schedules
fprintf('Running Oil-Water fully implicit model...\n');
tic;
[wsTFI, statesTFI] = simulateScheduleAD(state0, modelTFI, scheduleOW);
toc

fprintf('Running Oil-Water-Polymer fully implicit model...\n');
tic;
[wsPFI, statesPFI] = simulateScheduleAD(state0, modelPFI, scheduleP);
toc


%% Run sequential schedules
fprintf('Running Oil-Water sequential model...\n');
tic;
[wsTSQ, statesTSQ] = simulateScheduleAD(state0, modelTSQ, scheduleOW);
toc

%%
fprintf('Running Oil-Water-Polymer sequential model...\n');
tic;
[wsPSQ, statesPSQ] = simulateScheduleAD(state0, modelPSQ, scheduleP);
toc

%% Plot the accumulated water and oil production
wsTFIvc = vertcat(wsTFI{:});
wsPFIvc = vertcat(wsPFI{:});
wsTSQvc = vertcat(wsTSQ{:});
wsPSQvc = vertcat(wsPSQ{:});

qWs  = -([ [wsTFIvc(:,3).qWs] ; ...
           [wsPFIvc(:,3).qWs] ; ...
           [wsTSQvc(:,3).qWs] ; ...
           [wsPSQvc(:,3).qWs] ] .');
qWs  = bsxfun(@times, qWs, scheduleOW.step.val);
qOs  = -([ [wsTFIvc(:,3).qOs] ; ...
           [wsPFIvc(:,3).qOs] ; ...
           [wsTSQvc(:,3).qOs] ; ...
           [wsPSQvc(:,3).qOs] ] .');
qOs  = bsxfun(@times, qOs, scheduleOW.step.val);
cumt = cumsum(scheduleOW.step.val);

nCases = size(qWs, 2);
colors = lines(nCases);
linest = {'-','--','-','--'};

fh = figure;
op = get(fh, 'OuterPosition');
set(fh, 'OuterPosition', op.*[1 1 2 1]); %[left bottom width height]
hold on;
for i=1:nCases
	plot(convertTo(cumt, year), convertTo(qWs(:,i), stb), ...
    	'Color', colors(i,:), 'LineStyle', linest{i});
    plot(convertTo(cumt, year), convertTo(qOs(:,i), stb), ...
    	'Color', colors(i,:), 'LineStyle', linest{i});
end
legend({'Water, TFI', 'Oil, TFI', 'Water, PFI', 'Oil, PFI', ...
    'Water, TSQ', 'Oil, TSQ', 'Water, PSQ', 'Oil, PSQ'}, ...
	'Location', 'NorthEastOutside');
ylabel('Stb'); xlabel('Years');


%% Run POLYMER fully implicit schedule
fprintf('Running Oil-Water-Polymer fully implicit model...\n');
tic;
[wsPFI, statesPFI] = simulateScheduleAD(state0, modelPFI, scheduleP);
toc


%% Run POLYMER sequential schedule
fprintf('Running Oil-Water-Polymer sequential model...\n');
tic;
[wsPSQ, statesPSQ] = simulateScheduleAD(state0, modelPSQ, scheduleP);
toc


%% Plot the accumulated water and oil production
wsTFIvc = vertcat(wsTFI{:});
wsPFIvc = vertcat(wsPFI{:});
wsTSQvc = vertcat(wsTSQ{:});
wsPSQvc = vertcat(wsPSQ{:});

qWs  = -([ [wsTFIvc(:,3).qWs] ; ...
           [wsPFIvc(:,3).qWs] ; ...
           [wsTSQvc(:,3).qWs] ; ...
           [wsPSQvc(:,3).qWs] ] .');
qWs  = bsxfun(@times, qWs, scheduleOW.step.val);
qOs  = -([ [wsTFIvc(:,3).qOs] ; ...
           [wsPFIvc(:,3).qOs] ; ...
           [wsTSQvc(:,3).qOs] ; ...
           [wsPSQvc(:,3).qOs] ] .');
qOs  = bsxfun(@times, qOs, scheduleOW.step.val);
cumt = cumsum(scheduleOW.step.val);

nCases = size(qWs, 2);
colors = lines(nCases);
linest = {'-','-','--','--'};

fh = figure;
op = get(fh, 'OuterPosition');
set(fh, 'OuterPosition', op.*[1 1 2 1]); %[left bottom width height]
hold on;
for i=1:nCases
	plot(convertTo(cumt, year), convertTo(qWs(:,i), stb), ...
    	'Color', colors(i,:), 'LineStyle', linest{i});
    plot(convertTo(cumt, year), convertTo(qOs(:,i), stb), ...
    	'Color', colors(i,:), 'LineStyle', linest{i});
end
legend({'Water, TFI', 'Oil, TFI', 'Water, PFI', 'Oil, PFI', ...
    'Water, TSQ', 'Oil, TSQ', 'Water, PSQ', 'Oil, PSQ'}, ...
	'Location', 'NorthEastOutside');
ylabel('Stb'); xlabel('Years');

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

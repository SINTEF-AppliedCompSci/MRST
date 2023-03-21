%% Sleipner2019VEModel
% This example demonstrates how to setup a hybrid VE grid for the 2019 Sleipner
% benchmark model and run CO2 storage a simulation. This model is freely
% available as part of the CO2DataShare project
% 
% The full model is a layered model made up of 9 relatively thick sandstone
% layers divided by thin shale layers. In the hybrid VE model we represent
% the reservoir as a series of stacked VE models, one for each sandstone layer.
% We represent supposed high permeability chimneys which intersect the reservoir
% as fully resolved (non-VE) cells. There is no vertical flow modelled between
% the sandstone layers except at the chimneys.
%
% See "Watson et al. 2021. Rapid Optimisation of the new Sleipner Benchmark
% model" SINTEF Proceedings 7 Editors: Nils A. Røkke (SINTEF) and Hanna K. Knuutila (NTNU)
% TCCS–11CO2 Capture, Transport and Storage. Trondheim 22nd–23rd June 2021
% Short Papers from the 11th International Trondheim CCS Conference


%% Add modules
mrstModule add sleipner hybrid-ve co2lab wellpaths ad-core ad-props ad-blackoil deckformat coarsegrid libgeometry matlab_bgl mrst-gui

%% Get full grid
% The full model grid is read in from a .grdecl file. Within the function 
% getMultiLayerSleipnerGrid we also remove the caprock cells from model.
% The full model (including the caprock cells) can also be returned if 
% desired but we will not use it in this example.

[G, rock, grdecl, ~, ~] = getMultiLayerSleipnerGrid();

%% Get layer info
% This is needed to indicate which faces are at the top of each layer and used
% to define the edge of the VE layers.
disp('Getting layer boundaries.')
[layerBoundaryFaces] = getSleipnerLayerBoundaryFaces(G);
layerMapG = getLayerMapG(G);

%% Find Chimney cells and faces
% This will be used later to indicate which cells should be in the feeder
% chimneys.
disp('Finding feeder cells.')
[feeders] = getSleipnerFeederOutlines2019();
[feederCells,feederFaces,feederCellsMap] = findFeederIndicesInGrid(G,feeders,layerMapG);

%% Set permeabilities
% We convert all permeabilites to sand permeabilites to avoid any unwanted
% affects when upscaling layers into VE layers. The only exception is the thick
% shale layer near the top of the reservoir which is preserved as a VE layer
% with shale porosity and permeability.
disp('Setting permeabilties.')
rock.perm = rock.perm(:,1);

shalek = min(rock.perm);
sandk = max(rock.perm);

shaleporo = min(rock.poro);
sandporo = max(rock.poro);

% Set sand properties for all cells
rock.perm(:) = sandk;
rock.poro(:) = sandporo;

% Set shale permeabilities for thick shale layer 
rock.perm(layerMapG==9) = shalek;
rock.poro(layerMapG==9) = shaleporo;

% Set sand properties for all cells in fine chimneys
rock.perm(feederCells) = sandk;
rock.poro(feederCells) = sandporo;



%% Set fluid
% We get the fluid for the full model grid. Here were use an EOS based on P and T
% where we define the temperature in each cell at the start of the simulation
% and do not allow it to change.
disp('Getting fluid.')
topReservoirTemp = 37;
injTemp = 48;

[fluid] = getSleipner2019FluidModel('n',[1 1], ...
                                        'topReservoirTemp',topReservoirTemp,'injTemp',injTemp,...
                                        'G',G,'useEOS',true,'wcell',1596191); % first well cell in full model.

%% Make model for full grid
disp('Making base model.')
model= TwoPhaseWaterGasModel(G, rock, fluid);


%% Set transmissibilities
% We set zero transmissibility across faces at the layer boundaries.
% An exception is made for those faces within the feeder chimneys.
disp('Setting transmissibilities.')
tempFeederT = model.operators.T_all(feederFaces);
model.operators.T_all(layerBoundaryFaces) = 0;
model.operators.T_all(feederFaces) = tempFeederT;
model.operators.T = model.operators.T_all(model.operators.internalConn);

%% Make state0
disp('Calculating state0.')
region = getInitializationRegionsBlackOil(model, 0);
state0 = initStateBlackOilAD(model, region);

%% Setup wells for full grid
% Well trajectories are given in the dataset and are used to define wells in
% the model. We set a temporary injection rate which is modified later so that
% values match injection values given in Santi et al. (2018).
disp('Setting up schedule.')
injRate = 1*1000*mega./year; % This is temporary and will be modified later
injDensity = 760; % kg/m3 This is temporary and will be modified later

rate = injRate/injDensity;
W = addSleipnerWellsTrajectory(model.G,model.rock,rate);

%% Setup temporary schedule for full grid
ramp = rampupTimesteps(6, 6);
timesteps = [ramp(1:end-1);1;1;1; 1; 1; 1; 1; 1; 1]*year;
tmpschedule = simpleSchedule(timesteps, 'W', W, 'bc', []);


%% Modify schedule based on injection rates from Santi 2018 pg 22.
schedule.step.val = repmat(year,15,1);
schedule.step.control = (1:1:15)';

injRate = [0.07,0.66,0.84,0.93,0.93,1.01,0.96,0.91,0.75,0.86,0.82,0.92,0.81,0.86,0.74].*1000.*mega./year;

a = zeros(15,1);
for i  = 1:15
    schedule.control(i) = tmpschedule.control;
    schedule.control(i).W.val = injRate(i)/injDensity;
    a(i) = injRate(i)/injDensity;
end


%% Make VE model
disp('Making VE model.')
isFine = feederCells;
modelVE = convertToMultiVEModel(model, isFine, 'sealingFaces', layerBoundaryFaces);

%% Make new fluid using upscaled grid
% Because the EOS requires a temperature value at each cell we must recreate
% the fluid object with the VE grid.
disp('Making VE fluid.')
[fluidVE] = getSleipner2019FluidModel('n',[1 1], ...
                                          'topReservoirTemp',topReservoirTemp,'injTemp',injTemp,...
                                          'G',modelVE.G,'useEOS',true,'wcell',70645); %

modelVE.fluid = fluidVE;


%% Upscale schedule
scheduleVE = upscaleSchedule(modelVE, schedule, 'wellUpscaleMethod', 'sum', 'bcUpscaleMethod','mean');

%% Upscale state
disp('Upscaling state0.')
state0VE = upscaleState(modelVE, model, state0);

%% Add pressure bc
% We use pside function which generates a pressure boundary condition structure
% at the specified side of the grid. We won't use this structure but it is a 
% convenient way to get the required faces from the full grid.
bctemp = pside([],model.G,'east',10);
bctemp = pside(bctemp,model.G,'west',10);
bctemp = pside(bctemp,model.G,'north',10);
bctemp = pside(bctemp,model.G,'south',10);

% We then need to take the faces from the fine grid and identify which faces
% they correspond to in the VE grid.
CG = modelVE.G;
G = CG.parent;    
connCoarse = rldecode(1:CG.faces.num, diff(CG.faces.connPos), 2) .';

isFaceBC = false(G.faces.num, 1);
isFaceBC(bctemp.face) = true;

coarseFaceNo = zeros(G.faces.num, 1);
coarseFaceNo(CG.faces.fconn) = connCoarse;

coarseFacesBC = unique(connCoarse(isFaceBC(CG.faces.fconn)));

[bdryFaces,bdryCells] = boundaryFaces(modelVE.G);


cn = zeros(numel(coarseFacesBC),1);
for i = 1:numel(coarseFacesBC)
    cn(i) = bdryCells(bdryFaces==coarseFacesBC(i));
end

% Once we have the faces in the coarse grid we can find the hydrostatic pressure
% value used in the initial state and use this as the value at the boundary.
bcVal = state0VE.pressure(cn);
bc = addBC([], coarseFacesBC, 'pressure', bcVal, 'sat', [0 0]);

for i  = 1:15
    scheduleVE.control(i).bc = bc;
end

%% Make Problem
% We will use the packed problem setup to run the simulation and store results.
disp('Packing VE Problem.')
tic

pVE = packSimulationProblem(state0VE, modelVE, scheduleVE, 'Sleipner2019VE');

%% Simulate
simulatePackedProblem(pVE);

%% Visualise results
% First convert the VE results into results on the fine grid
[~, states, ~] = getPackedSimulatorOutput(pVE);
states_finescale = convertMultiVEStates(pVE.SimulatorSetup.model, states);
% Plot results interactively using the plotToolbar function
plotToolbar(pVE.SimulatorSetup.model.G.parent, states_finescale)

%% Plot CO2 saturation at the top of each layer at the end of the simulation

figure()
for i = 1:9
    subplot(3,3,i)
    if i == 9
        lix = 10;
    else
        lix = i;        
    end

    cells = (layerMapG==lix);
    plotCellData(pVE.SimulatorSetup.model.G.parent,states_finescale{end}.s(:,2),...
                 cells,'EdgeColor','none');
    hold on
    title(strcat("L",num2str(i)))
end

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.
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

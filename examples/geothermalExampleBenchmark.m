%% Benchmark with TOUGH2
% In this example we simulate injection of pure water at 50C in a
% reservoir of brine at 10C. The same simulation was done with TOUGH2
% The results are compared to validate the implementation of equations in
% the geothermal module. 
%
% Setup: 2D grid without gravity -> equivalent to a top view 
% Pressure applied on the left side (injection), on two faces only
% Reservoir with constant pressure, temperature and NaCl mass fraction
% conditions (100 bars, 10C, 0.1)
% 
% Pressure, temperature and NaCl mass fraction are evaluated at two cells
% in the model: near the injection and in the model centre.

%% Add necessary MRST modules
mrstModule add ad-props ad-core ad-blackoil compositional geothermal mrst-gui

%% Set up the grid
physdim  = [100 200];              % domain size in x,y directions 
celldim  = [50 100];               % number of cells in x, y dims 
G = cartGrid(celldim, physdim);    % Cartesian Grid 
G = computeGeometry(G);            % grid geometry and connectivity 

%% Find injection faces and monitoring cell indices
% monitoring point near injection
Minj = find(G.cells.centroids(:,1)==1 & G.cells.centroids(:,2)==101);   
% Monitoring point in the center
Mcen = find(G.cells.centroids(:,1)==51 & G.cells.centroids(:,2)==101);   
% Injection boundary faces
Finj1 = find(G.faces.centroids(:,1)==0 & G.faces.centroids(:,2)==99);
Finj2 = find(G.faces.centroids(:,1)==0 & G.faces.centroids(:,2)==101);
faces = [Finj1; Finj2];

%% Plot the grid with the monitoring cells
plotting_grid = true; % set to true to plot the grid
if plotting_grid
figure('Position', [0,0,450,700])
plotGrid(G, 'faceColor', 'none', 'edgeAlpha', 0.2)
hold on 
plotGrid(G, Minj, 'faceColor', 'green')
plotGrid(G, Mcen, 'faceColor', 'blue')
axis equal 
axis tight
title('Grid and location of monitoring cells')
txt1 = 'Injection';
txt2 = 'Centre';
text(5, 95, txt1,'FontSize',12)
text(55, 105, txt2,'FontSize',12)
end 

%% Set fluid structure properties
rhoWS = 1000;
% Define fluid structure
fluid = initSimpleADIFluid('mu'    , 1.0e-3, ...
                           'rho'   , rhoWS , ...
                           'phases', 'W'   );
fluid = addThermalFluidProps(fluid           , ... % Original fluid
                             'Cp'     , 4.2e3, ... % Specific heat capacity
                             'lambdaF', 0.6  , ... % Thermal conductivity
                             'useEOS' , true );    % Use equation of state

%% Make rock structure
perm = 1e-14;
poro = 0.1;
% define rock structure
rock = makeRock(G, perm, poro);
% add thermal props
rock = addThermalRockProps(rock           , ... % Original rock
                           'CpR'    , 1000, ... % Specific heat capacity
                           'lambdaR', 2   , ... % Thermal conductivity
                           'rhoR'   , 2700, ... % Rock density
                           'tau'    , 1   );    % Tortuosity

%% Define the model
% Provide data for one-phase two-component model with H2O and NaCl
compFluid = CompositionalBrineFluid(          ...
    {'H2O'             , 'NaCl'            }, ... % Names
    [18.015281*gram/mol, 58.442800*gram/mol], ... % Molar masses
    [0                 , 1e-6              ]);    % Molecular diffusivities
% Make model
model = GeothermalModel(G, rock, fluid, compFluid, 'extraStateOutput', true);

%% Set up initial state
state0   = initResSol(G, 100*barsa, 1);       % pressure and saturation
state0.T = ones(G.cells.num,1).*283.15;       % temperature
X = repmat([0.9, 0.1], G.cells.num, 1);       % Initial mass fractions
state0.components = model.getMoleFraction(X); % Convert to mole fractions

%% Set Boundary conditions and schedule
bc  = pside([], G, 'right', 100*barsa, 'sat', 1);        % Right (p)
bc  = addBC(bc, faces, 'pressure', 250*barsa ,'sat', 1); % Left (p)
Tbc = [repmat(283.15, 100, 1); repmat(323.15, 2, 1)];
bc  = addThermalBCProps(bc, 'T', Tbc);                   % Temperature
Xbc = [repmat(0.1, 100, 1); zeros(2, 1)];
bc.components = model.getMoleFraction([1-Xbc, Xbc]);     % Mole fractions

% Define the same timesteps as for the simulation with Tough2
fn   = fullfile(getDatasetPath('geothermal'), 'tough2-results.mat');
data = load(fn); tough2 = data.tough2;
time      = tough2(:,2);
timesteps = diff([0;time]);
schedule  = simpleSchedule(timesteps, 'bc', bc);

%% Run simulation
[~, states] = simulateScheduleAD(state0, model, schedule);

%% Post processing
% Get data at injection monitoring point
model = model.validateModel();
pinj = cellfun(@(st) st.pressure(Minj), states);
Tinj = cellfun(@(st) st.T(Minj), states);
Xinj = cell2mat(cellfun(@(st) model.getMassFraction(st.components(Minj,:)), states, 'UniformOutput', false));
Xinj = Xinj(:,2);
% ... and centre monitoring point
pcen = cellfun(@(st) st.pressure(Mcen), states);
Tcen = cellfun(@(st) st.T(Mcen), states);
Xcen = cell2mat(cellfun(@(st) model.getMassFraction(st.components(Mcen,:)), states, 'UniformOutput', false));
Xcen = Xcen(:,2);

%% Plot benchmark tough vs. MRST
tough2Ix   = [4,9,7,12,5,10];
tough2Data = tough2(:, tough2Ix);
tough2Data(:,1:2) = tough2Data(:,1:2)./barsa;
tough2Data(:,3:4) = tough2Data(:,3:4).*0.1;

mrstData = [pinj/barsa , pcen/barsa , ...
            Xinj       , Xcen       , ...
            Tinj-273.15, Tcen-273.15];

plotTime = cumsum(schedule.step.val)/day;

yl = {'Pressure [bar]', 'NaCl mass fraction', 'Temperature [C]'};
tl = {'Injection monitoring point', 'Centre monitoring point'};

figure('Position', [0,0,1200,800])
for i = 1:numel(tough2Ix)
    subplot(3,2,i), hold on
    plot(plotTime, tough2Data(:,i), '-' , 'LineWidth', 2);
    plot(plotTime, mrstData(:,i)  , '--', 'LineWidth', 5);
    hold off, box on
    ylabel(yl{ceil(i/2)});
    xlabel('Time [day]');
    title(tl{mod(i-1,2)+1})
    if i == 1
        legend({'TOUGH2', 'MRST'},'Location','northeast')
    end
end

%% Copyright Notice
%
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
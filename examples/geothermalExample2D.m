%% Simple 2D case with gravity and variable fluid properties 
% Simulate mass and heat transfer in a reservoir. The model considers
% injection at the left boundary and outflow at the right boundary
% (pressure set at reservoir condition). The top and bottom layer have
% no-flow boundary conditions. The initial reservoir is set at 100 bars and
% 10 C. Pure water is injected at 50 C.

%% Add necessary MRST modules
mrstModule add ad-props ad-core ad-blackoil compositional geothermal mrst-gui

%% Set plot options
setAxProps = @(ax) set(ax, 'View'              , [-140,20]    , ...
                           'PlotBoxAspectRatio', [2,2,1]      , ...
                           'Projection'        , 'Perspective', ...
                           'Box'               , 'on'         , ...
                           'XLimSpec'          , 'tight'      , ...
                           'YLimSpec'          , 'tight'      , ...
                           'ZLimSpec'          , 'tight'      );

%% Set up the grid
physdim  = [100 50];               % domain size in x, y directions
celldim  = [50  50];               % number of cells in x, y directions
G = cartGrid(celldim, physdim);    % Cartesian grid
G = computeGeometry(G);            % grid geometry and connectivity

%% Make fluid structure properties
% Define fluid structure
fluid = initSimpleADIFluid('mu', 1, 'rho', 1, 'phases', 'W');
% Add thermal properties. The module supports using an EOS from Spivey et
% al (2004) with p/T-dependent density and viscosity. For convenience, this
% can be added to the fluid by setting the optional argument 'useEOS' to
% true in 'addThermalFluidProps'.
fluid = addThermalFluidProps(fluid, 'Cp'     , 4.2e3, ...
                                    'lambdaF', 0.6  , ...
                                    'useEOS' , true );

%% Inspect fluid model
% We plot the fluid density and viscosity as a function of pressure and
% temperature. We also indicate the normal operational range by a black
% triangle
K0   = 273.15*Kelvin;
pMin = 1e6*Pascal;      % Minimum pressure
pMax = 200e6*Pascal;    % Maximum pressure
TMin = K0;              % Minimum temperature 
TMax = K0 + 275*Kelvin; % Maximum temperature 
n    = 100;
% Get pressure/temperature grid
[p,T] = ndgrid(linspace(pMin, pMax, n), linspace(TMin, TMax, n));
% "Normal" operational range
p1    = linspace(pMin, 500*barsa, n)';
T1    = ones(n, 1)*(K0 + 100*Kelvin);
p2    = ones(n, 1)*500*barsa;
T2    = linspace(TMin, K0 + 100*Kelvin, n)';
% Plot density
rhoW = reshape(fluid.rhoW(p(:),T(:)), n, n);
figure(), surf(p/barsa, T-K0, rhoW, 'EdgeColor', 'none'); setAxProps(gca);
xlabel('Pressure (bar)'), ylabel('Temperature (C)'), zlabel('Density, (kg/m^3)');
hold on
rhoW1 = fluid.rhoW(p1, T1);
plot3(p1/barsa, T1-K0, rhoW1, 'k');
rhoW2 = fluid.rhoW(p2, T2);
plot3(p2/barsa, T2-K0, rhoW2, 'k');
hold off
% Plot viscosity
muW = reshape(fluid.muW(p(:),T(:)), n, n);
figure(), surf(p/barsa, T-K0, muW/(centi*poise), 'EdgeColor', 'none'); setAxProps(gca);
xlabel('Pressure (bar)'), ylabel('Temperature (C)'), zlabel('Viscosity, (cp)');
hold on
muW1 = fluid.muW(p1, T1);
plot3(p1/barsa, T1-K0, muW1/(centi*poise), 'k');
muW2 = fluid.muW(p2, T2);
plot3(p2/barsa, T2-K0, muW2/(centi*poise), 'k');
hold off

%% Make rock structure
perm = 1e-14;
poro = 0.1;  
% Define rock structure
rock = makeRock(G, perm, poro);   
% Add thermal props
rock = addThermalRockProps(rock, 'lambdaR', 2, 'rhoR', 2700, 'CpR', 1000);

%% Inititate state transient flow
state0   = initResSol(G, 100*barsa, 1);       % Pressure and saturation 
state0.T = ones(G.cells.num,1).*(273.15+10);  % Temperature

%% Define the problem - GeothermalWaterModel for pure water problem
% Define the gravity vector for a 2D x,y case using MRST convention
gravity reset on
gravity([0 -9.81]); % gravity is downward 
model = GeothermalModel(G, rock, fluid, 'extraStateOutput', true);
% wModel  =  GeothermalWaterModel(G, rock, fluid, 'extraStateOutput', true);
% Define limits of temperature validity for EOS
model.minimumTemperature = TMin;
model.maximumTemperature = TMax;
model.maximumPressure    = pMax;

%% Set up boundary conditions
pRes = 100*barsa;
pInj = pRes + 100*barsa;
Tinj = K0+50;
Tres = K0+10;
% Flow BCs
bc  = [];
bc  = pside(bc, G, 'left', pInj, 'sat', 1); nf   = numel(bc.face);
bc  = pside(bc, G, 'right', pRes, 'sat', 1);
Tbc = repmat(Tres, numel(bc.face), 1); Tbc(1:nf) = Tinj;
bc  = addThermalBCProps(bc, 'T', Tbc);

%% Set up schedule
timesteps = rampupTimesteps(1*year, 30*day, 8); 
schedule  = simpleSchedule(timesteps, 'bc', bc);

%% Run simulation
[~, states] = simulateScheduleAD(state0, model, schedule);

%% Visualization
figure(), plotToolbar(G, states); axis tight equal, colormap(hot)

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
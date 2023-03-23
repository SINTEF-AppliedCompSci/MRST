close all;
clear;

%% add the necessary modeules 
mrstModule add ad-core ad-props ad-blackoil geochemistry mrst-gui
mrstVerbose off

%% Define the domain 
% here we look at a 2D cartesian grid  with 100 cells

sideLength = 10;
G = cartGrid([sideLength, sideLength, 1], [sideLength, sideLength, 1]);
G = computeGeometry(G);
nc = G.cells.num;

% you can view the domain using plotGrid
plotGrid(G), view(3), axis tight
outInd = [G.faces.num-G.cells.num-sideLength+1; G.faces.num-G.cells.num-sideLength^2+sideLength];

% here we plot the faces of a zero dirchlet boundary condition in pressure
plotFaces(G, outInd, 'r');

%% Define the rock
% porosity and permeability are the minimum number of parameters for rock

rock.perm = 1*darcy*ones(G.cells.num, 1);
rock.poro = 0.5*ones(G.cells.num, 1);

%% Define the fluid
pRef = 0*barsa;
fluid = initSimpleADIFluid('phases', 'W', 'mu', 1*centi*poise, 'rho', ...
                           1000*kilogram/meter^3, 'c', 1e-10, 'cR', 4e-10, ...
                           'pRef', pRef);

%% Define the chemistry
% here we look at a simple tracer chemical system 

elements = {'O', 'H', 'Na*','Cl*'};

species = {'H+*', 'OH-', 'Na+', 'H2O*', 'NaCl','Cl-'};

reactions ={'H2O  = H+  + OH- ',      10^-14*mol/litre, ... 
            'NaCl = Na+ + Cl-',       10^1*mol/litre};


% instantiate chemical model
chemsys = ChemicalSystem(elements, species, reactions);

chemmodel = ChemicalModel(chemsys);

%% solve the chemistry for the initial and injected compositions
% we will be injected a low slainity fluid into a high salinity aquifer. we
% solve the system so that we can give the total element concentrations
% within the domain at the initial time, and define the boundary conditions

% initial chemistry
Nai = 1e-1;
Cli = Nai;
Hi = 1e-7;
H2Oi = 1;
initial = [Nai Cli Hi H2Oi]*mol/litre;

% we must repeat the initial chemistry for each cell of the system. This
% can be done by passing a vector to initState, or repmatting the state
% variable produced
[initChemState, initreport]= chemmodel.initState(repmat(initial, nc,1), 'charge', 'Cl');

% the initial state must also contain a pressure field 
initChemState.pressure = pRef*ones(nc,1);

% injected chemistry
Naf = 1e-3;
Clf = Naf;
Hf = 1e-10;
H2Of = 1;
injected = [Naf Clf Hf H2Of]*mol/litre;
[injChemState, injreport]= chemmodel.initState(injected, 'charge', 'Cl');


%% Define the transport model
model = ChemicalTransportModel(G, rock, fluid, chemmodel);

%% Define the boundary conditions
% here we have two source cells in the southwest and northeast corners.
% there are two dirchlet conditions in the northwest and southeast corners
% allowing for outflow. Defining only neumann boundary conditions will
% underconstain the pressure solution such that a constant of integration
% is unknown.

% grab the volume of each cell in the domain
pv = poreVolume(G,rock);

% inject one pore volume into the cells per day.
src                	= [];
src               	= addSource(src, [1; nc], pv(1:2)/day, 'sat', 1);

% define chemistry of the injected fluid
src.elements        = injChemState.elements(end,:);
src.logElements     = injChemState.logElements(end,:);

% specifiy dirchlet conditions for the composition
bc                  = [];
outInd = [G.faces.num-G.cells.num-sideLength+1; G.faces.num-G.cells.num-sideLength^2+sideLength];
bc                  = addBC(bc, outInd, 'pressure', [0; 0]*barsa, 'sat', 1);
bc.elements         = initChemState.elements(end,:);        % will not used if outflow
bc.logElements      = initChemState.logElements(end,:);  % will not used if outflow


%% Define the time stepping schedule
% here we define the time stepping schedule. it is recommended to ramp up
% the size of the time steps gradually

schedule.step.val = [0.01*day*ones(5, 1); 0.1*day*ones(5,1); 1*day*ones(5, 1); 5*day*ones(10, 1)];
schedule.step.control = ones(numel(schedule.step.val), 1);
schedule.control = struct('bc', bc, 'src', src, 'W', []);

%% Run the simulation

[~, states, scheduleReport] = simulateScheduleAD(initChemState, model, schedule);

%% visualize the simulation
plotToolbar(G, states, 'field', 'pressure','startplayback', true)

%% Copyright notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2017 SINTEF ICT, Applied Mathematics.
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


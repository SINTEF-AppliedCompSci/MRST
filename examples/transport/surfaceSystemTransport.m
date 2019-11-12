close all;
clear;

mrstModule add ad-core ad-props ad-blackoil geochemistry

mrstVerbose off

%% Define the grid
% here we define a 1D domain with 10 cells, each a meter in length

G = cartGrid([10, 1, 1], [10, 1, 1]);
G = computeGeometry(G);
nc = G.cells.num;

%% Define the rock

rock.perm = 1*darcy*ones(G.cells.num, 1);
rock.poro = ones(G.cells.num, 1);

%% Define the fluid
pRef = 0*barsa;
fluid = initSimpleADIFluid('phases', 'W', 'mu', 1*centi*poise, 'rho', ...
                           1000*kilogram/meter^3, 'c', 1e-10, 'cR', 4e-10, ...
                           'pRef', pRef);

%% Define the chemistry
% the chemical system include the speciation of NaCl and H2O, as well as a
% reactive silica surface which can have a negative or neutral surface

elements = {'O', 'H', 'Na*','Cl*'};

species = {'H+*', 'OH-', 'Na+', 'H2O*', '>SiO-', '>SiOH', 'NaCl','Cl-'};

reactions ={'H2O  = H+  + OH- ',      10^-14*mol/litre, ... 
            'NaCl = Na+ + Cl-',       10^1*mol/litre,...
            '>SiOH = >SiO- + H+',     10^-8*mol/litre};

geometry = [2*site/(nano*meter)^2 50e-3*meter^2/(gram) 5e3*gram/litre];
sioInfo = {geometry, 'tlm', [1 1e3], '>SiO-', [-1 0 0], '>SiOH', [0 0 0]};

surfaces = {'>SiO', sioInfo};

% instantiate chemical model object
chemsys = ChemicalSystem(elements, species, reactions, 'surf', surfaces);

% print the chemical model to the screen
chemsys.printChemicalSystem;

chemmodel = ChemicalModel(chemsys);

%% solve for the initial state
% the column is originally saturated with a basic solution, an acidic
% solution is injected on the left boundary

% initial chemistry
Nai = 1e-3;
Cli = Nai;
Hi = 1e-9;
H2Oi = 1;

inputConstraints = [Nai Cli Hi H2Oi]*mol/litre;
[initchemstate, initreport]= chemmodel.initState(repmat(inputConstraints, nc,1), 'charge', 'Cl');

initState = initchemstate;
initState.pressure = pRef*ones(nc,1);

% injected chemistry
Naf = 1e-1;
Clf = Nai;
Hf = 1e-9;
H2Of = 1;

inputConstraints = [Naf Clf Hf H2Of]*mol/litre;
[injchemstate, injreport] = chemmodel.initState(inputConstraints, 'charge', 'Cl');

%% Define the transport model
model = ChemicalTransportModel(G, rock, fluid, chemmodel);

%% Define the boundary conditions
% here we specify a dirchlet boundary condition for pressure on the right
% hand side, and a constant flux on the left. The total element
% concentration at the boundaries must be provided.

% use model.fluidMat to pull the fluid concentrations from the injected
% state
injfluidpart  = injchemstate.species*model.fluidMat';
initfluidpart = initchemstate.species(end,:)*model.fluidMat';

pv = poreVolume(G,rock);

% define source term at cell 1 for inflow boundary condition
src = [];
src = addSource(src, 1, pv(1)/day, 'sat', 1);

% give the fluid concentration at the inlet
src.elements = injfluidpart;
src.logElements = log(injfluidpart);

% give dirchlet boundary condition at outlet
bc             = [];
bc             = pside(bc, G, 'east', 0*barsa, 'sat', 1);
bc.elements    = initfluidpart;        % (will not be used if outflow)
bc.logElements = log(initfluidpart);  % (will not be used if outflow)


%% Define the schedule
% it is recommened to ramp up the time stepping

% ten time steps of 0.01 days followed by 100 steps of 1 day
schedule.step.val     = [0.01*day*ones(10, 1); 0.1*day*ones(10, 1); 1*day*ones(10, 1);];
schedule.step.control = ones(numel(schedule.step.val), 1);
schedule.control      = struct('bc', bc, 'src', src, 'W', []);


%% Run the simulation

[~, states, scheduleReport] = simulateScheduleAD(initState, model, schedule);

%% visualize the simulation
mrstModule add mrst-gui
plotToolbar(G, states, 'field', 'species:1','startplayback', true, 'plot1D', true)

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


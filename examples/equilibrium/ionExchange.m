clear;
close all;

%% add adi and geochemistry modules
mrstModule add geochemistry ad-core
mrstVerbose on
%% instantiate the chemical model
% here we look at the competetive sorption of sodium and protons for an
% ionexhcnage surface. ChemicalModel will verify that surfaces that are
% designated as ion exchange surfaces do not contain species with a charge

elements = {'O', 'H', 'Na*', 'Cl*'};

species = {'H+*', 'OH-', 'H2O*', '>YH', '>YNa', 'Na+', 'Cl-','NaOH'};

reactions = {'H2O  = H+  + OH- ',           10^-14*mol/litre, ...
             '>YH + Na+ = >YNa + H+',       10^1,...
             'NaOH = Na+ + OH-',            10^10*mol/litre};

surfInfo = {'>Y', {[1*site/(nano*meter)^2 1*meter^2/gram 1*gram/litre], 'ie'}};

chemsys = ChemicalSystem(elements, species, reactions, 'surfaces', surfInfo);

% print the chemical system
chemsys.printChemicalSystem;

% Setup model
chemmodel = ChemicalModel(chemsys);

%% specify inputs
% here we vary pH
n = 100;

H2O = ones(n,1);
H   = logspace(-3, -11,n)';
Na  = 1e-2*ones(n,1);
Cl  = Na;

inputs = [Na, Cl, H, H2O]*mol/litre;

%% solve the chemical system
% here we specify that charge balance should be enforced by giving the
% name/value pair. The total amount of sodium will be varied in order for
% the charge balance equation to be satisfied

[state, report, chemmodel] = chemmodel.initState(inputs, 'chargeBalance', 'Na');

%% process the data
state = changeUnits(state, {'elements','species'}, mol/litre );

[state, chemmodel] = chemmodel.updateActivities(state);

pH = -log10(getProp(chemmodel, state, 'aH+'));

%% plot the results
figure; hold on; box on;
plot(pH, state.species, 'linewidth', 2)

xlabel('pH')
ylabel('concentration [mol/L]');
set(gca, 'yscale', 'log');
legend(chemsys.speciesNames)
xlim([3 11])

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
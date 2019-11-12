clear;
close all;

%% load ad and geochemistry modules
mrstModule add geochemistry ad-core
mrstVerbose on

%% generate chemical system object
% here we look at a simple example of aqueous speciation  of sodium,
% protons, and chlorine

elements = {'O', 'H', 'Na*', 'Cl*'};

species = {'H+*', 'OH-', 'Na+', 'Cl-', 'NaCl', 'H2O*'};

reactions = {'H2O  = H+  + OH- ', 10^-14*mol/litre, ...
             'NaCl = Na+ + Cl-',  10^1*mol/litre};

chemsys = ChemicalSystem(elements, species, reactions);

chemsys.printChemicalSystem;

%% define the input 
n = 100;

Na  = ones(n,1)*milli*mol/litre;
Cl  = ones(n,1)*milli*mol/litre;
H2O = ones(n,1)*mol/litre;
H   = logspace(-4, -10, n)'*mol/litre;

inputs = [Na, Cl, H, H2O];


%% solve the chemical system
chemmodel = ChemicalModel(chemsys);
[state, report, chem] = chemmodel.initState(inputs, 'chargeBalance', 'Na');

%% process the data
state = changeUnits(state, {'elements','species'}, mol/litre );

[state, chem] = chemmodel.updateActivities(state);

pH = -log10(getProp(chem, state, 'aH+'));

%% plot the results
figure; hold on; box on;
plot(pH, state.species, 'linewidth',2)

xlabel('pH')
ylabel('concentration [mol/L]');
set(gca, 'yscale', 'log');
legend(chemsys.speciesNames)

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
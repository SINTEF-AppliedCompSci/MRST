clear;
close all;

%% load adi and geochemistry module
mrstModule add ad-core geochemistry
mrstVerbose on
%% generate chemical system
% here we look at boron sorption onto the silica surface. we are keeping
% Na, Cl, B, H+ and H2O concentrations as inputs

% define elements
elements = {'O', 'H', 'Na*', 'Cl*', 'B*'};

species = {'H+*', 'OH-', 'Na+', 'Cl-', 'NaCl', 'H2O*',...
             'H3BO3', 'H4BO4-',...
             '>SO-' '>SOH', '>SOH2+','>SOH3BO3-'};

reactions ={'H2O  = H+  + OH- ',              10^-14*mol/litre, ...
            'NaCl = Na+ + Cl-',               10^1*mol/litre,...
            'H3BO3 + H2O = H4BO4- + H+',      10^-9.2,...
            '>SOH + H+ = >SOH2+',             10^8.31/(mol/litre),...
            '>SOH = >SO- + H+',               10^-11.82*mol/litre,...
            '>SOH + H3BO3 = >SOH3BO3- + H+'    10^-8.2};

% define surface information. Here we describe the surface with the
% constant capacitance model with a capacitance of 1.06 after the values in
% Goldberg et al 2000
geometry = [2.31*site/(nano*meter)^2 59*meter^2/gram 200*gram/litre];
sInfo = {geometry, 'ccm', 1.06};
surfaces ={ '>SO', sInfo};

% instantiate chemical model
chemsys = ChemicalSystem(elements, species, reactions, 'surf', surfaces);

% print the chemical system
chemsys.printChemicalSystem;

% Setup model
chemmodel = ChemicalModel(chemsys);

%% define input parameters
% here we vary the pH of the system
n = 200;

Na  = 1e-1.*ones(n,1);
Cl  = 1e-1.*ones(n,1);
B   = 4.63e-4.*ones(n,1);
H2O = ones(n,1);
H   = logspace(-7,-11, n)';

userInput = [Na Cl B H H2O]*mol/litre;

%% solve the chemical system
[state, report, model] = chemmodel.initState(userInput, 'charge', 'Cl');


%% calculate auxillary information
% surface charge, potential, aqueous and surface concentrations can be
% calculated with tools built into ChemicalModel

[state, chemmodel] = chemmodel.updateActivities(state);
[state, chemmodel] = chemmodel.updateChargeBalance(state);
[state, chemmodel] = chemmodel.updateSurfacePotentials(state);
[state, chemmodel] = chemmodel.updateAqueousConcentrations(state);
[state, chemmodel] = chemmodel.updateSurfaceConcentrations(state);


state = changeUnits(state, {'species','activities','elements','surfaceConcentrations'}, mol/litre);
pH = -log10(getProp(chemmodel, state, 'aH+'));

%% data from Goldberg
% data is from Goldberg et al. Soil Sci. Soc. Am. J. 64:1356-1363 (2000)
GB =   [7.145895	0.2151590
        7.259432	0.2151301
        7.277191	0.2502701
        7.399525	0.2578790
        7.662131	0.3296291
        7.977074	0.3937258
        8.029960	0.4517771
        8.380156	0.5540653
        8.669038	0.6349768
        9.097366	0.6807084
        9.176045	0.6898564
        9.604094	0.7019716
        9.813089	0.6285734
        10.29323	0.6040028
        10.65009	0.4572219
        10.81594	0.4464835];


%% plot the results
figure; hold on; box on;
plot(pH, state.chargeBalance)
xlabel('pH')
ylabel('charge [% of total ion concentration]');

% aqueous components
figure; hold on; box on;
plot(pH, chemmodel.getProp(state, 'H3BO3'));
plot(pH,  chemmodel.getProp(state, 'H4BO4-'));
legend('H_3BO_3','H_4BO_4^-');

xlabel('pH');
ylabel('concentration [mol/L]');

% surface components
figure; hold on; box on;
plot(pH, chemmodel.getProp(state, '>SOH3BO3-')/200*1e6, '-k')
plot(GB(:,1), GB(:,2), 'ok');
legend('mrst', 'data')
   
xlabel('pH')
ylabel('adsorbed B [\mu mol/g]');

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


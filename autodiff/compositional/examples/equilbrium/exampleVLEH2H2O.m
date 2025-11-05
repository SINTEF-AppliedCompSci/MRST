%% H2-H2O Vapor-Liquid Equilibrium Calculations
% This MRST example calculates the vapor-liquid equilibrium (VLE) for a hydrogen-water
% mixture under varying conditions using a thermodynamic equation of state.
% Two EoS models are used: Peng-Robinson and Soreide-Whitson.
% Then we compare the Solubility of H2 in water and NaCl brine obtained
% with both EoS models
% Author: [Stéphanie Delage Santacreu]
% Date: [16/09/2025]
% Organization: [Université de Pau et des Pays de l'Adour, E2S UPPA, CNRS, LFCR, UMR5150, Pau, France]
% ---------------------------------------------------------------------------

clear; clc;

% Add necessary MRST modules
mrstModule add h2-biochem compositional ad-blackoil ad-core ad-props mrst-gui

%% Define Composition Mixture
% Create a compositional mixture with water and hydrogen components
compFluid = TableCompositionalMixture({'Water', 'Hydrogen'}, {'H2O', 'H2'});
disp(compFluid);

%% Initialize Thermodynamic Model
% Choose an equation of state model for the calculations
eosNamesw = 'sw'; % Soreide-Whitson (SW) model
eosModelsw = EquationOfStateModel([], compFluid, eosNamesw);

eosNamepr = 'pr'; % Peng Robinson (PR) model
eosModelpr = EquationOfStateModel([], compFluid, eosNamepr);

%% Define Test Case Parameters
% Set the test case for different pressures, temperatures, and salinity levels
z0 = [0.8, 0.2]; % Initial composition
caseTest = 1; % Choose the test case here

switch caseTest
    case 1 %sources: Solubility of H2 in water and NaCl brine under subsurface
        % storage conditions (https://hal.science/hal-04623907v1, 2023)
        eosModelsw.msalt=0;
        pres=[100.01, 150.01, 200.01, 200.01, 101.11, 101.31, 130.01,...
            165.01, 199.91, 200.11, 100.01, 100.01, 100.01, 175.11]*barsa;
        Temp=[298.20,298.05,298.15,298.15,323.55,323.50,323.85,323.55,...
            323.30,323.35,373.85,373.80,373.85,373.65]*Kelvin;
        xliqH2Exp=[0.00135994,0.00199679,0.00264397,0.00263320,0.00125091,...
            0.00123925,0.00159253,0.00201117,0.00244186,0.00245479,...
            0.00142471,0.00140332,0.00140893,0.00243347];

    case 2
        eosModelsw.msalt=1;
        Temp=[298.20,298.30,298.15,298.30,323.20,323.40,323.35,323.20,...
            323.30,323.30,323.20,373.25,373.40,373.10,373.15,373.45,373.15,373.00];
        pres=[100.71,150.01,150.01,200.01,100.31,100.61,101.01,150.06,...
            175.01,199.91,200.01,100.11,126.01,150.36,150.46,150.71,175.51,200.46]*barsa;
        xliqH2Exp=[0.00107012,0.00159298,0.00161244,0.00213575,0.00102099,...
            0.00102078,0.00102426,0.00151377,0.00176728,0.00202590,...
            0.00204487,0.00119671,0.00148492,0.00175595,0.00179573,...
            0.00177350,0.00204000,0.00234604];
    case 3
        eosModelsw.msalt=2;
        Temp=[298.15,298.05,298.05,323.20,323.40,323.35,323.40,373.05,373.20,373.40];
        pres=[100.01,150.01,200.01,100.01,150.01,150.01,200.01,100.01,150.01,200.51]*barsa;
        xliqH2Exp=[0.00088640,0.00132235,0.00171848,0.00088260,0.00131242,...
            0.00128912,0.00172402, 0.00099379, 0.00151866, 0.00205031 ];
    case 4
        eosModelsw.msalt=4;
        Temp=[298.20, 298.20, 298.15, 323.30, 323.40, 323.40, 373.25, 373.35, 373.15];
        pres=[100.01, 150.01, 200.01, 100.01, 150.01, 200.01, 100.01, 150.01, 200.01]*barsa;
        xliqH2Exp=[0.00059422, 0.00093595, 0.00121838, 0.00061736, ...
            0.00095752, 0.00129237, 0.00077991, 0.00114509, 0.00157469];
end


%% Perform Flash Calculations
% Determine the liquid-phase hydrogen fraction (xliqH2) for each condition
nc = numel(pres);
namecp = eosModelsw.getComponentNames();
indH2=find(strcmp(namecp,'H2'));
indH2O= find(strcmp(namecp,'H2O'));

[Lsw, xsw, ~] = standaloneFlash(pres, Temp, z0, eosModelsw);
[Lpr, xpr, ~] = standaloneFlash(pres, Temp, z0, eosModelpr);
xliqH2sw=xsw(:,indH2);
xliqH2pr=xpr(:,indH2);



%% calculate and display the errors models vs experiments
errorSWexp=abs(xliqH2sw-xliqH2Exp')./xliqH2Exp';
errorPRexp=abs(xliqH2pr-xliqH2Exp')./xliqH2Exp';

fprintf('Errormax SW-Experiment: %12.8f, Errormax PR-Experiment: %12.8f\n', max(errorSWexp), max(errorPRexp));
fprintf('Errormin SW-Experiment: %12.8f, Errormin PR-Experiment: %12.8f\n', min(errorSWexp), min(errorPRexp));
fprintf('Errormean SW-Experiment: %12.8f, Errormean PR-Experiment: %12.8f\n', mean(errorSWexp), mean(errorPRexp));

%% plot the results
plotEosPRSW(pres,xliqH2Exp,xliqH2sw,xliqH2pr)

%% Copyright notice

% <html>
% <p><font size="-1">
% Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.
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
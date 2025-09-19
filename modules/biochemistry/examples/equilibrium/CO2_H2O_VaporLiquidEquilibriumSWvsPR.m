%% C02-H2O Vapor-Liquid Equilibrium Calculations
% This script calculates the vapor-liquid equilibrium (VLE) for a CO2-water
% mixture under varying conditions using a thermodynamic equation of state.
% Two EoS models are used: Peng-Robinson and Soreide-Whitson.
% Then we compare the Solubility of CO2 in water and NaCl brine obtained 
% with both EoS models 
% Author: [Stéphanie Delage Santacreu]
% Date: [16/09/2025]
% Organization: [Université de Pau et des Pays de l'Adour, E2S UPPA, CNRS, LFCR, UMR5150, Pau, France]

clear; clc;

% Add necessary MRST modules
mrstModule add biochemistry compositional ad-blackoil ad-core ad-props mrst-gui

%% Define Composition Mixture
% Create a compositional mixture with water and carbon dioxide components
compFluid = TableCompositionalMixture({'Water', 'CarbonDioxide'}, {'H2O', 'CO2'});
disp(compFluid);
%% Initialize Thermodynamic Model
% Choose an equation of state model for the calculations
eosNamesw = 'sw'; % Soreide-Whitson (SW) model
eosModelsw = SoreideWhitsonEquationOfStateModel([], compFluid, eosNamesw);

eosNamepr = 'pr'; % Peng Robinson (PR) model
eosModelpr = SoreideWhitsonEquationOfStateModel([], compFluid, eosNamepr);

%% Define Test Case Parameters
% Set the test case for different pressures, temperatures, and salinity levels
z0 = [0.8, 0.2]; % Initial composition
caseTest = 2; % Choose the test case here

switch caseTest %Source: Thermodynamic study of the CO2-H2O system, S.Chabab (https://hal.science/hal-02310963v1,2019)
     case 1       
        eosModel.msalt=1.13;
        Temp=[323.02, 322.97, 323.03, 323.04, 372.33, 372.31, 372.29, 372.29, 372.25]*Kelvin;
        pres=[53.450, 75.550, 100.350, 145.080, 31.148, 60.500, 108.840, 151.920, 191.980]*barsa;
        xliqCO2Exp=[0.01030, 0.01290, 0.01510, 0.01700,0.00390, 0.00750, 0.01130, 0.01360, 0.01570];
            
       
     case 2
        eosModel.msalt=1;
        Temp=[373.38, 373.37, 373.41]*Kelvin;
        pres=[16.983, 32.527, 68.182]*barsa;
        xliqCO2Exp=[0.00237, 0.00426, 0.00833];
        
    case 3
        eosModel.msalt=3.01;
        Temp=[342.82, 342.81, 342.82, 372.39, 372.42, 372.41 , 372.43, 372.45, 372.45]*Kelvin;
        pres=[30.391, 72.559, 100.910, 25.556, 71.417, 100.517, 152.433, 199.597, 229.817]*barsa;
        xliqCO2Exp=[0.00441, 0.00880, 0.01057, 0.00292, 0.00707, 0.00878, 0.01141, 0.01258, 0.01337];

    case 4
        eosModel.msalt=0;
        Temp=[323.2, 323.2, 323.2, 323.2, 323.2, 323.2, 333.2, 333.2, 333.2,...
            333.2, 333.2, 333.2, 353.1, 353.1, 353.1, 353.1, 353.1, 353.1]*Kelvin;
        pres=[40.5, 60.6, 80.8, 100.9, 121, 141.1, 40.5, 60.6, 80.8, 100.9,...
            121, 141.1, 40.5, 60.6, 80.8, 100.9, 121, 131]*barsa;
        xliqCO2Exp=[0.0109, 0.0161, 0.019, 0.0205, 0.0214, 0.0217, 0.0096,...
            0.0138, 0.0166, 0.0186, 0.0201, 0.0208, 0.008, 0.0114, 0.014, ...
            0.016, 0.0176, 0.0184];

end


%% Perform Flash Calculations
% Determine the liquid-phase hydrogen fraction (xliqH2) for each condition
nc = numel(pres);
namecp = eosModelsw.getComponentNames();
indCO2=find(strcmp(namecp,'CO2'));
indH2O= find(strcmp(namecp,'H2O'));

[Lsw, xsw, ~] = standaloneFlash(pres, Temp, z0, eosModelsw); 
[Lpr, xpr, ~] = standaloneFlash(pres, Temp, z0, eosModelpr); 
xliqCO2sw=xsw(:,indCO2);
xliqCO2pr=xpr(:,indCO2);



%% calculate and display the errors models vs experiments
error_swexp=abs(xliqCO2sw-xliqCO2Exp')./xliqCO2Exp';
error_prexp=abs(xliqCO2pr-xliqCO2Exp')./xliqCO2Exp';

fprintf('Errormax SW-Experiment: %12.8f, Errormax PR-Experiment: %12.8f\n', max(error_swexp), max(error_prexp));
fprintf('Errormin SW-Experiment: %12.8f, Errormin PR-Experiment: %12.8f\n', min(error_swexp), min(error_prexp));
fprintf('Errormean SW-Experiment: %12.8f, Errormean PR-Experiment: %12.8f\n', mean(error_swexp), mean(error_prexp));

%% plot the results
plot_EosPRSW(pres,xliqCO2Exp,xliqCO2sw,xliqCO2pr)

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
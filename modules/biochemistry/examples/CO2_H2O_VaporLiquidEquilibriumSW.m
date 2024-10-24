%% CO2-H2O Vapor-Liquid Equilibrium Calculations
% This script calculates the vapor-liquid equilibrium (VLE) for a CarbonDioxide-water
% mixture under varying conditions using a thermodynamic equation of state.

clear; clc;

% Add necessary MRST modules
mrstModule add biochemistry compositional ad-blackoil ad-core ad-props mrst-gui

%% Define Composition Mixture
% Create a compositional mixture with water and hydrogen components
compFluid = TableCompositionalMixture({'Water', 'CarbonDioxide'}, {'H2O', 'CO2'});
disp(compFluid);

%% Initialize Thermodynamic Model
% Choose an equation of state model for the calculations
eosName = 'sw'; % Soreide-Whitson (SW) model
eosModel = SoreideWhitsonEquationOfStateModel([], compFluid, eosName);

%% Define Test Case Parameters
% Set the test case for different pressures, temperatures, and salinity levels
z = [0.8, 0.2]; % Initial composition
patm = 1e5; % Atmospheric pressure in Pa
caseTest =4; % Choose the test case here

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
% Determine the liquid-phase hydrogen fraction (xliqCO2) for each condition
nc = numel(pres);
namecp = eosModel.getComponentNames();
indCO2=find(strcmp(namecp,'CO2'));
indH2O= find(strcmp(namecp,'H2O'));
[L, x, ~] = standaloneFlash(pres, Temp, z, eosModel); 
xliqCO2=x(:,indCO2);

%% Write Results to File
% Save the results (temperature, pressure, hydrogen mole fraction) to a file
outputFileName = sprintf('CO2solubility_case%d_msalt%d_%s.dat', caseTest, eosModel.msalt, eosName);
fileID = fopen(outputFileName, 'wt');

fprintf(fileID, 'Temperature (K)  Pressure (bar)  CO2_molar_fraction  CO2_molar_fraction_exp\n');
for i = 1:nc
    fprintf(fileID, '%12.2f  %12.4f  %12.8f  %12.8f\n', Temp(i), pres(i) / patm, xliqCO2(i), xliqCO2Exp(i));
end
fclose(fileID);

disp(['Results written to file: ', outputFileName]);



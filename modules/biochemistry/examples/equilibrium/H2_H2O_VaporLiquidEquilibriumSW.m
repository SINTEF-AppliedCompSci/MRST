%% H2-H2O Vapor-Liquid Equilibrium Calculations
% This script calculates the vapor-liquid equilibrium (VLE) for a hydrogen-water
% mixture under varying conditions using a thermodynamic equation of state.

clear; clc;

% Add necessary MRST modules
mrstModule add biochemistry compositional ad-blackoil ad-core ad-props mrst-gui

%% Define Composition Mixture
% Create a compositional mixture with water and hydrogen components
compFluid = TableCompositionalMixture({'Water', 'Hydrogen'}, {'H2O', 'H2'});
disp(compFluid);

%% Initialize Thermodynamic Model
% Choose an equation of state model for the calculations
eosName = 'sw'; % Soreide-Whitson (SW) model
eosModel = SoreideWhitsonEquationOfStateModel([], compFluid, eosName);

%% Define Test Case Parameters
% Set the test case for different pressures, temperatures, and salinity levels
z = [0.8, 0.2]; % Initial composition
patm = 1e5; % Atmospheric pressure in Pa
caseTest = 4; % Choose the test case here

switch caseTest
    case 1 %sources: Solubility of H2 in water and NaCl brine under subsurface 
        % storage conditions (https://hal.science/hal-04623907v1, 2023)
        eosModel.msalt=0;
        pres=[100.01, 150.01, 200.01, 200.01, 101.11, 101.31, 130.01,...
            165.01, 199.91, 200.11, 100.01, 100.01, 100.01, 175.11]*barsa;     
        Temp=[298.20,298.05,298.15,298.15,323.55,323.50,323.85,323.55,...
            323.30,323.35,373.85,373.80,373.85,373.65]*Kelvin;
        xliqH2Exp=[0.00135994,0.00199679,0.00264397,0.00263320,0.00125091,...
            0.00123925,0.00159253,0.00201117,0.00244186,0.00245479,...
            0.00142471,0.00140332,0.00140893,0.00243347];

    case 2
        eosModel.msalt=1;
        Temp=[298.20,298.30,298.15,298.30,323.20,323.40,323.35,323.20,...
            323.30,323.30,323.20,373.25,373.40,373.10,373.15,373.45,373.15,373.00];
        pres=[100.71,150.01,150.01,200.01,100.31,100.61,101.01,150.06,...
            175.01,199.91,200.01,100.11,126.01,150.36,150.46,150.71,175.51,200.46]*barsa;
        xliqH2Exp=[0.00107012,0.00159298,0.00161244,0.00213575,0.00102099,...
            0.00102078,0.00102426,0.00151377,0.00176728,0.00202590,...
            0.00204487,0.00119671,0.00148492,0.00175595,0.00179573,...
            0.00177350,0.00204000,0.00234604];
    case 3
        eosModel.msalt=2;
        Temp=[298.15,298.05,298.05,323.20,323.40,323.35,323.40,373.05,373.20,373.40];
        pres=[100.01,150.01,200.01,100.01,150.01,150.01,200.01,100.01,150.01,200.51]*barsa;
        xliqH2Exp=[0.00088640,0.00132235,0.00171848,0.00088260,0.00131242,...
            0.00128912,0.00172402, 0.00099379, 0.00151866, 0.00205031 ];
    case 4
        eosModel.msalt=4;
        Temp=[298.20, 298.20, 298.15, 323.30, 323.40, 323.40, 373.25, 373.35, 373.15];
        pres=[100.01, 150.01, 200.01, 100.01, 150.01, 200.01, 100.01, 150.01, 200.01]*barsa;
        xliqH2Exp=[0.00059422, 0.00093595, 0.00121838, 0.00061736, ...
            0.00095752, 0.00129237, 0.00077991, 0.00114509, 0.00157469];
end

%% Perform Flash Calculations
% Determine the liquid-phase hydrogen fraction (xliqH2) for each condition
nc = numel(pres);
namecp = eosModel.getComponentNames();
indH2=find(strcmp(namecp,'H2'));
indH2O= find(strcmp(namecp,'H2O'));
[L, x, ~] = standaloneFlash(pres, Temp, z, eosModel); 
xliqH2=x(:,indH2);

%% Write Results to File
% Save the results (temperature, pressure, hydrogen mole fraction) to a file
outputFileName = sprintf('H2solubility_case%d_msalt%d_%s.dat', caseTest, eosModel.msalt, eosName);
fileID = fopen(outputFileName, 'wt');

fprintf(fileID, 'Temperature (K)  Pressure (bar)  H2_molar_fraction  H2_molar_fraction_exp\n');
for i = 1:nc
    fprintf(fileID, '%12.2f  %12.4f  %12.8f  %12.8f\n', Temp(i), pres(i) / patm, xliqH2(i), xliqH2Exp(i));
end
fclose(fileID);

disp(['Results written to file: ', outputFileName]);



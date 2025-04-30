function getFluidH2BrineProps(tab_h2o, tab_h2, tab_sol, varargin)
% getFluidH2BrineProps - Write fluid property tables based on gas dissolution state.
%
% This function computes the thermodynamic properties for a fluid system
% consisting of hydrogen (H2) and brine (water with dissolved salts). It
% utilizes input data tables for pure water, pure hydrogen, and solubility
% data to determine properties such as density, viscosity, vapor-liquid
% ratios, and mole fractions. The calculated results are written to output
% text files for subsequent analysis or simulation use.
%
% Syntax:
%   getFluidH2BrineProps(tab_h2o, tab_h2, tab_sol, varargin)
%
% Input Parameters:
%   tab_h2o (table):
%       A structured dataset containing thermodynamic properties of pure water.
%       It should include parameters such as temperature, density, viscosity,
%       and saturation properties essential for calculations involving the water phase.
%
%   tab_h2 (table):
%       A structured dataset containing thermodynamic properties of pure hydrogen.
%       This dataset must encompass relevant parameters like density and viscosity
%       at varying pressures and temperatures, which are critical for understanding
%       hydrogen behavior in the mixture.
%
%   tab_sol (array or table):
%       A structured dataset representing the solubility of hydrogen in brine.
%       It should include information on how hydrogen concentration changes with
%       pressure and temperature, along with any relevant properties for calculating
%       the phase behavior of the system.
%
%   varargin (optional):
%       A variable-length input argument list that allows for additional options
%       and configurations. This can include parameters such as file output directory
%       or specific flags to control the function's behavior, enabling more flexible usage.
%
% Output:
%   The function does not return outputs directly. Instead, it writes calculated
%   properties to text files, including:
%   - PVTO.txt: Contains pressure, solution gas ratio, water density, and viscosity
%     for the brine phase.
%   - PVTG.txt: Contains pressure, vapor-liquid ratio, gas density, and viscosity
%     for the hydrogen phase.
% SEE ALSO:
%   generateSolubilityTable, generateComponentTable.

%{
Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.

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

% Import necessary modules
mrstModule add deckformat
%% Define default options
% Default options
opts = struct('rs', true, ...
    'rv', true, ...
    'plot', true, ...
    'nusat', 10, ...
    'outputPath', '', ...
    'units', 'metric', ...
    'reCompute', true, ...
    'PVTGFile', 'PVTG_UHS.txt', ...
    'PVDOFile', 'PVDO_UHS.txt', ...
    'PVTOFile', 'PVTO_UHS.txt');

% Merge user-defined options
opts = merge_options(opts, varargin{:});

% Ensure the output directory exists
if isempty(opts.outputPath)
    opts.outputPath = fullfile(mrstOutputDirectory(), 'UHS_PVT', 'PVT_H2STORAGE');
end

if ~exist(opts.outputPath, 'dir')
    mkdir(opts.outputPath);
end
%% Configuration parameters
do_plot = opts.plot;
dissolve_gas = opts.rs;
vaporize_water = opts.rv;

%% Extract and align pressures from the tables based on temperature
T = tab_h2.("# temperature [°C]")(1);
tab_sol = filterTableByTemperature(tab_sol, T);
tab_h2o = filterTableByTemperature(tab_h2o, T);
tab_h2 = filterTableByTemperature(tab_h2, T);

%% Extract pressure values for consistency check
p_sol = tab_sol.("pressure [Pa]");
p_h2o = tab_h2o.("pressure [Pa]");
p_h2 = tab_h2.("pressure [Pa]");

% Define pressure range and temperature
pressureRange = [min(p_sol), max(p_sol)];
fprintf('Tabulating PVT data for pressure range: [%.2f, %.2f] Pa at temperature: %.2f °C\n', ...
    pressureRange(1), pressureRange(2), T);

% Assert that the filtered tables are non-empty
assert(~isempty(tab_sol) && ~isempty(tab_h2o) && ~isempty(tab_h2), 'Filtered tables are empty.');
assert(height(tab_sol) == height(tab_h2o) && height(tab_sol) == height(tab_h2), 'Mismatch in table sizes.');

% Check if pressures are consistent within a tolerance of 0.001 Pa
tolerance = 1e-3;
assert(all(abs(p_sol - p_h2o)./barsa() < tolerance), 'Error: Pressure mismatch detected between H2O and solution.');
assert(all(abs(p_sol - p_h2)./barsa() < tolerance), 'Error: Pressure mismatch detected between H2 and solution.');

%% Convert pressure to bar
p = p_sol ./ barsa;

%% Handle dissolved gas and vaporized water phases
[x_H2, y_h2o] = computeMoleFractions(tab_sol, dissolve_gas, vaporize_water);
[X_H2, Y_h2o] = convertMoleToMassFraction(x_H2, y_h2o);
%% Plot phase properties if plotting is enabled
if do_plot
    figure;
    subplot(1, 2, 1);
    plot(p_sol, x_H2, 'LineWidth', 2);
    xlabel('Pressure (Pa)');
    ylabel('x_{H2} (Liquid Mole Fraction)');
    title('x_{H2} in Liquid');

    subplot(1, 2, 2);
    plot(p_sol, y_h2o, 'LineWidth', 2);
    xlabel('Pressure (Pa)');
    ylabel('y_{H2O} (Vapor Mole Fraction)');
    title('y_{H2O} in Vapor');
end

%% Retrieve densities and viscosities
[pure_water_density, pure_gas_density, water_viscosity, gas_viscosity] = extractFluidProperties(tab_h2o, tab_h2);

%% Plot densities if enabled
if do_plot
    figure; hold on;
    saturated_water_density = WaterDensity(T, pure_water_density, X_H2);

    subplot(1, 2, 1);
    plot(p, pure_water_density, 'LineWidth', 2, 'DisplayName', 'Pure Water');
    hold on;
    plot(p, saturated_water_density, 'LineWidth', 2, 'DisplayName', 'H2 Saturated Water');
    xlabel('Pressure (Pa)');
    ylabel('Density (kg/m^3)');
    title('Water Density');
    legend;
    hold off;

    subplot(1, 2, 2);
    plot(p, pure_gas_density, 'LineWidth', 2);
    xlabel('Pressure (Pa)');
    ylabel('Density (kg/m^3)');
    title('Gas Density');

    hold off;
end

%% Surface densities for calculation
rhoGS = 0.0851; % Density of H2 at surface (kg/m^3)
rhoOS = 9.98207150467e+02; % Density of H2O at surface (kg/m^3)

%% Calculate saturation properties Rs and Rv
if dissolve_gas
    R_s = rhoOS .* X_H2 ./ (rhoGS .* (1 - X_H2));
end

if vaporize_water
    R_v = rhoGS .* Y_h2o ./ (rhoOS .* (1 - Y_h2o));
end

if do_plot
    figure;
    % Plot R_s if dissolve_gas is true
    if dissolve_gas
        subplot(2, 1, 1); % 2 rows, 1 column, first plot
        plot(p_sol, R_s, 'LineWidth', 2);
        title('Saturated R_s');
        xlabel('Pressure (Pa)');
        ylabel('R_s (kg/m³)');
        grid on;
    end

    % Plot R_v if vaporize_water is true
    if vaporize_water
        subplot(2, 1, 2); % 2 rows, 1 column, second plot
        plot(p_sol, R_v, 'LineWidth', 2);
        title('Saturated R_v');
        xlabel('Pressure (Pa)');
        ylabel('R_v (kg/m³)');
        grid on;
    end
end

%% Compute shrinkage factors for liquid and gas
bO = shrinkageFactor(WaterDensity(T, pure_water_density, X_H2), Y_h2o, rhoOS, rhoGS);
bG = shrinkageFactor(pure_gas_density, Y_h2o, rhoGS, rhoOS);

if do_plot
    figure; hold on;

    % Plot b_g (Gas Formation Volume Factor)
    subplot(1, 2, 1);
    plot(p, bG, 'LineWidth', 2);
    title('Saturated b_g');
    xlabel('Pressure (Pa)');
    ylabel('b_g (m³/Sm³)');
    grid on;

    % Plot b_o (Oil Formation Volume Factor)
    subplot(1, 2, 2);
    plot(p, bO, 'LineWidth', 2);
    title('Saturated b_o');
    xlabel('Pressure (Pa)');
    ylabel('b_o (m³/Sm³)');
    grid on;
end


%% Write fluid property tables based on dissolution state
if dissolve_gas
    writePVTO(p_sol, R_s, pure_water_density, water_viscosity, T, X_H2, Y_h2o, rhoOS, rhoGS, opts);
else
    writeIMMISC(p_sol, pure_water_density ./ rhoOS, water_viscosity, 'PVDO', opts);
end
if vaporize_water
    writePVTG(p_sol, R_v, pure_gas_density, gas_viscosity, Y_h2o, rhoGS, rhoOS, opts);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function writeIMMISC(p, b, mu, title, opts)
% Function to write immiscible properties to a text file
% Inputs:
%   p     - Pressure values (vector)
%   b     - Formation volume factor (vector)
%   mu    - Viscosity values (vector)
%   title - Title/name for the output file (string)
%   opts  - Options structure containing directory path and units

% Construct the output file path
filePath = fullfile(opts.outputPath, opts.PVDOFile);
fileName = opts.PVDOFile;
if exist(filePath, 'file')&&(~opts.reCompute)
    fprintf('File %s already exists. Loading data...\n', fileName); 
    return; % Exit the function as the table already exists
end

fn = fopen(filePath, 'w');

% Check if the file was opened successfully
if fn == -1
    error('Cannot open file: %s', filePath);
end

% Get unit conversion factors based on specified units
u = unitConversionFactors('si', opts.units);

% Determine the formation volume factor unit based on title
if strcmpi(title, 'pvdg')
    fvf_u = u.gasvol_r / u.gasvol_s;  % Gas volume factor
else
    fvf_u = u.liqvol_r / u.liqvol_s;  % Liquid volume factor
end

% Write header for the output file
fprintf(fn, '-- PRES    FVF      VISC\n');
fprintf(fn, '%s\n', title);  % Write the title

% Loop through each pressure value and write data to file
for i = 1:numel(p)
    fprintf(fn, '%-4.10g %-4.10g %-4.10g\n', ...
        p(i) * u.press, ...        % Pressure in proper units
        fvf_u / b(i), ...         % Formation volume factor
        u.viscosity * mu(i));     % Viscosity in proper units
end

fprintf(fn, '/\n');  % End of the data block
fclose(fn);          % Close the file

% Display confirmation message
fprintf('%s written!\n', filePath);
end

function writePVTO(p, Rs, pure_water_density, water_viscosity, T, X_h2_sat, Y_h2o_sat, rhoOS, rhoGS, opts)
% Function to write PVTO data to a file.
% Inputs:
% p: Pressure values (vector)
% Rs: Solution gas-to-oil ratio (vector)
% pure_water_density: Density of pure water (vector)
% water_viscosity: Viscosity of water (vector)
% T: Temperature (scalar)
% X_h2_sat: Saturation of hydrogen (vector)
% Y_h2o_sat: Saturation of water (vector)
% rhoOS: Oil phase density (scalar)
% rhoGS: Gas phase density (scalar)
% n: Number of saturation pressure points
% opts: Options structure containing directory path and units

% Prepare file path and open the file for writing
filePath = fullfile(opts.outputPath, opts.PVTOFile); % Use the filename from opts
fileName = opts.PVTOFile;
if exist(filePath, 'file')&&(~opts.reCompute)
    fprintf('File %s already exists. Loading data...\n', fileName);
    return;
end

fn = fopen(filePath, 'w');
nusat = opts.nusat;


% Check if the file was opened successfully
if fn == -1
    error('Cannot open file: %s', filePath);
end

% Write header for the file
fprintf(fn, '-- RS    PRES    FVF      VISC\n');
fprintf(fn, 'PVTO\n');

% Preprocessing for minimum conditions
X_h2_min = 1.0e-8; % Minimum hydrogen saturation
Rs_min = rhoOS * X_h2_min / rhoGS; % Minimum Rs calculation
p_zero = interp1(X_h2_sat, p, X_h2_min, 'linear', 'extrap');
assert(p_zero > 0, 'Pressure at minimum hydrogen saturation must be positive.');

% Extend pressure and saturation values
rhoW_min = interp1(p, pure_water_density, p_zero, 'linear', 'extrap');
p = [p_zero; p];                     % Extend pressure vector
Rs = [Rs_min; Rs];                   % Extend Rs vector
X_h2_sat = [X_h2_min; X_h2_sat];     % Extend X_h2 saturation
Y_h2o_sat = [Y_h2o_sat(1); Y_h2o_sat]; % Extend Y_h2o saturation
pure_water_density = [rhoW_min; pure_water_density]; % Extend water density
water_viscosity = [water_viscosity(1); water_viscosity]; % Extend water viscosity

% Compute shrinkage factor
bO_sat = shrinkageFactor(WaterDensity(T, pure_water_density, X_h2_sat), X_h2_sat, rhoOS, rhoGS);

% Loop through each Rs and write data
M = numel(Rs);
u = unitConversionFactors('si', opts.units); % Get unit conversion factors
uk = u.gasvol_s / u.liqvol_s;   % Conversion factor from gas volume to liquid volume
ub = u.liqvol_r / u.liqvol_s;   % Conversion factor for liquid volume in reservoir conditions

for i = 1:M
    mu = water_viscosity(i); % Current water viscosity
    % Write Rs, pressure, FVF, and viscosity to file
    fprintf(fn, '%-4.10g %-4.10e %-4.10g %-4.10g\n', ...
        Rs(i) * uk, p(i) * u.press, ub * (1 / bO_sat(i)), u.viscosity * mu);

    % Compute pressure values for unsaturated conditions
    X_h2 = X_h2_sat(i);
    p_usat = linspace(p(i), 1.2 * max(p), nusat + 1); % Unsaturated pressure values
    p_usat = p_usat(2:end); % Exclude first element (equal to p(i))

    % Write data for each unsaturated pressure point
    for j = 1:nusat
        rho_w_j = interp1(p, pure_water_density, p_usat(j), 'linear', 'extrap'); % Interpolate water density
        rho_w_dissolved = WaterDensity(T, rho_w_j, X_h2); % Calculate density with dissolved hydrogen
        bO_usat = shrinkageFactor(rho_w_dissolved, X_h2, rhoOS, rhoGS); % Calculate shrinkage factor
        rho_new = bO_usat * (rhoOS + Rs(i) * rhoGS); % New density calculation
        assert(abs(rho_w_dissolved - rho_new) < 1e-3, 'Density mismatch exceeds tolerance.');

        % Write unsaturated pressure, FVF, and viscosity values to file
        fprintf(fn, '           %-4.10g %-4.10g %-4.10g', ...
            p_usat(j) * u.press, ub * (1 / bO_usat), u.viscosity * mu);
        if j == nusat
            fprintf(fn, ' /\n'); % End of line for last entry
        else
            fprintf(fn, '\n'); % New line for subsequent entries
        end
    end
end

% Close the file
fclose(fn);
fprintf('%s written!\n', filePath); % Confirmation message
end


function writePVTG(p, Rv, pure_gas_density, gas_viscosity, Y_h2o, rhoGS, rhoOS, opts)
% Function to write PVTG data to a file
% Inputs:
% p: Pressure values (vector)
% Rv: Vaporization ratios (vector)
% pure_gas_density: Density of the pure gas (scalar)
% gas_viscosity: Viscosity of the gas (vector)
% Y_h2o: Water content (scalar)
% rhoGS: Gas phase density (scalar)
% rhoOS: Oil phase density (scalar)
% opts: Options structure containing directory path and units

% Prepare file path and open the file for writing
filePath = fullfile(opts.outputPath, opts.PVTGFile);
fileName = opts.PVTGFile;
if exist(filePath, 'file')&&(~opts.reCompute)
    fprintf('File %s already exists. Loading data...\n', fileName); 
    return; % Exit the function as the table already exists
end
fn = fopen(filePath, 'w');

% Check if the file was opened successfully
if fn == -1
    error('Cannot open file: %s', filePath);
end

% Write header for the file
fprintf(fn, '-- PRES       RV       FVF       VISC\n');
fprintf(fn, 'PVTG\n');

% Calculate shrinkage factor once before using it
bG = shrinkageFactor(pure_gas_density, Y_h2o, rhoGS, rhoOS);
% Get conversion factors
u = unitConversionFactors('si', opts.units);
uk = u.press;           % Pressure conversion factor
urv = u.liqvol_r / u.gasvol_s;  % Liquid volume to gas volume conversion factor
ub = u.gasvol_r / u.gasvol_s;   % Gas volume conversion factor
umu = u.viscosity;     % Viscosity conversion factor

% Loop through each Rv and write data
M = numel(Rv);
for i = 1:M
    mu = gas_viscosity(i); % Gas viscosity for current index
    % Now you can index into bG
    fprintf(fn, '%-10.5g %-10.12g %-10.8g %-10.12g\n', ...
        p(i) * uk, Rv(i) * urv, ub / bG(i), mu * umu);
    fprintf(fn, '             %-10.5g %-10.8g %-10.12g\n', 0.0, ub / bG(i), mu * umu);


    fprintf(fn, ' /\n'); % End of data block
end

% Close the file
fclose(fn);
fprintf('%s written!\n', filePath); % Confirmation message
end



function rho = HydrogenDensityInWater(T)
% Computes the density of hydrogen in water based on temperature.
% Inputs:
% T: Temperature in Celsius (scalar)
% Outputs:
% rho: Density of hydrogen in water (kg/m^3)

% Calculate the partial volume of H2 in water at temperature T
V = PartialVolumeOfH2InWater(T);
M = 2.016e-3; % Molar mass of H2 (kg/mol)

% Calculate the density of H2
rho = M ./ V; % Density (kg/m^3)
end

function V = PartialVolumeOfH2InWater(T)
% Computes the partial volume of hydrogen in water based on temperature.
% Inputs:
% T: Temperature in Celsius (scalar)
% Outputs:
% V: Partial volume of hydrogen in water (m^3/mol)

% Coefficients for the partial volume calculation
a = 51.1904; % Coefficient a
b = -2.08062e-1; % Coefficient b
c = 3.4427e-4; % Coefficient c

% Convert temperature to Kelvin
T_K = T + 273.15;

% Calculate partial volume of H2 in water
V = (a + b .* T_K + c .* T_K.^2) ./ 1e6; % Convert from cm^3/mol to m^3/mol
end

function rhoW = WaterDensity(T, pure_water_density, X_h2)
% Computes the density of a water-hydrogen mixture based on temperature and composition.
% Inputs:
% T: Temperature in Celsius (scalar)
% pure_water_density: Density of pure water (kg/m^3)
% X_h2: Mole fraction of hydrogen (scalar)
% Outputs:
% rhoW: Density of the mixture (kg/m^3)

% Calculate the mole fraction of water
X_h2o = 1 - X_h2;

% Compute the volume contributions from water and hydrogen
vol = X_h2o ./ pure_water_density + X_h2 ./ HydrogenDensityInWater(T);

% Calculate the density of the mixture
rhoW = 1 ./ vol; % Density (kg/m^3)
end

function bO = shrinkageFactor(rhoO, X_h2, rhoOS, rhoGS)
% Computes the shrinkage factor for a mixture of oil and gas based on density and composition.
% Inputs:
% rhoO: Density of the oil phase (kg/m^3)
% X_h2: Mole fraction of hydrogen (scalar)
% rhoOS: Density of oil-saturated phase (kg/m^3)
% rhoGS: Density of gas-saturated phase (kg/m^3)
% Outputs:
% bO: Shrinkage factor (dimensionless)

% Calculate the ratio of oil saturation to gas saturation
Rs = rhoOS .* X_h2./ (rhoGS.*(1-X_h2));

% Calculate the shrinkage factor
bO = rhoO ./ (rhoOS + Rs * rhoGS); % Dimensionless factor
end

function mf = moleFractionToMassFraction(x, molar_mass_self, molar_mass_other)
% Converts mole fraction to mass fraction.
% Calculate the mass fraction
mf = (x * molar_mass_self) ./ (x .* molar_mass_self + (1 - x) .* molar_mass_other); % Mass fraction
end

function [X_H2, Y_h2o] = convertMoleToMassFraction(x_H2, y_h2o)
    molar_mass_h2 = 2.016e-3;
    molar_mass_h2o = 18.01528e-3;
    X_H2 = moleFractionToMassFraction(x_H2, molar_mass_h2, molar_mass_h2o);
    Y_h2o = moleFractionToMassFraction(y_h2o, molar_mass_h2o, molar_mass_h2);
end

function [x_H2, y_h2o] = computeMoleFractions(tab_sol, dissolve_gas, vaporize_water)
x_H2 = 0;
y_h2o = 0;
if dissolve_gas
    x_H2 = max(tab_sol.x_H2 .* dissolve_gas, 1e-8);
end
if vaporize_water
    y_h2o = tab_sol.y_H2O .* vaporize_water;
end
end

function table_filtered = filterTableByTemperature(table, T)
    table_filtered = table(abs(table.("# temperature [°C]") - T) < eps, :);
end

function [rho_w, rho_g, mu_w, mu_g] = extractFluidProperties(tab_h2o, tab_h2)
    rho_w = tab_h2o.('density [kg/m3]');
    rho_g = tab_h2.('density [kg/m3]');
    mu_w = tab_h2o.('viscosity [Pa.s]');
    mu_g = tab_h2.('viscosity [Pa.s]');
end

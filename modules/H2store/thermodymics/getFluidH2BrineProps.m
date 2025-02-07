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

% Import necessary modules
mrstModule add deckformat
%% Define default options
% Default options
    opts = struct('rs', true, ...
                   'rv', true, ...
                   'plot', false, ...
                   'nusat', 10, ...
                   'dir', '', ...
                   'units', 'metric', ...
                   'PVTGFile', 'PVTG.txt', ...
                   'PVDOFile', 'PVDO.txt', ...
                   'PVTOFile', 'PVTO.txt');

    % Merge user-defined options
    opts = merge_options(opts, varargin{:});
%% Set directory for output files; create if it does not exist
if ~isempty(opts.dir) && ~isfolder(opts.dir)
    mkdir(opts.dir);
end

%% Assert that the directory exists now
assert(isfolder(opts.dir), 'The specified directory could not be created or does not exist.');

%% Configuration parameters
do_plot = opts.plot;
dissolve_gas = opts.rs;
vaporize_water = opts.rv;

% Extract temperature
T = tab_h2.x_Temperature__C_(1); % Assuming reservoir temperature is the first

%% Extract and align pressures from the tables based on temperature
tab_sol = tab_sol(abs(tab_sol.x_Temperature__C_ - T) < eps, :);
tab_h2o = tab_h2o(abs(tab_h2o.x_Temperature__C_ - T) < eps, :);
tab_h2 = tab_h2(abs(tab_h2.x_Temperature__C_ - T) < eps, :);


%% Extract pressure values for consistency check
p_sol = tab_sol.phasePressure_Pa_;
p_h2o = tab_h2o.pressure_Pa_;
p_h2 = tab_h2.pressure_Pa_;

% Define pressure range and temperature
pressureRange = [min(p_sol), max(p_sol)];
fprintf('Tabulating PVT data for pressure range: [%.2f, %.2f] Pa at temperature: %.2f °C\n', ...
    pressureRange(1), pressureRange(2), T);

% Assert that the filtered tables are non-empty
assert(~isempty(tab_sol), 'No matching temperature found in tab_sol for the specified T.');
assert(~isempty(tab_h2o), 'No matching temperature found in tab_h2o for the specified T.');
assert(~isempty(tab_h2), 'No matching temperature found in tab_h2 for the specified T.');
assert(height(tab_sol) == height(tab_h2o) && height(tab_sol) == height(tab_h2), ...
    'Filtered tables do not have the same x-size (number of pressure points); check pressure range.');

% Check if pressures are consistent within a tolerance of 0.001 Pa
tolerance = 0.001;
assert(all(abs(p_sol - p_h2o) < tolerance), 'Mismatch in pressures between H2O and solution.');
assert(all(abs(p_sol - p_h2) < tolerance), 'Mismatch in pressures between H2 and solution.');

%% Convert pressure to bar
p = p_sol ./ barsa;

%% Constants for H2 and H2O molar masses
molar_mass_h2 = 2.016e-3; % kg/mol
molar_mass_h2o = 18.01528e-3; % kg/mol

%% Handle dissolved gas and vaporized water phases
x_H2 = max(tab_sol.x_H2___ .* dissolve_gas, 1e-8);
y_h2o = tab_sol.y_H2O___ .* vaporize_water;
X_H2 = mole_fraction_to_mass_fraction(x_H2, molar_mass_h2, molar_mass_h2o);
Y_h2o = mole_fraction_to_mass_fraction(y_h2o, molar_mass_h2o, molar_mass_h2);

%% Plot phase properties if plotting is enabled
if do_plot
    figure(1); clf;
    subplot(1, 2, 1); plot(p_sol, x_H2); title('x_H2 in liquid');
    subplot(1, 2, 2); plot(p_sol, y_h2o); title('y_H2O in vapor');
end

%% Retrieve densities and viscosities
pure_water_density = tab_h2o.density_kg_m3_;
pure_gas_density = tab_h2.density_kg_m3_;
water_viscosity = tab_h2o.viscosity_Pa_s_;
gas_viscosity = tab_h2.viscosity_Pa_s_;

%% Compute saturated water density
saturated_water_density = water_density(T, pure_water_density, X_H2);

% Plot densities if enabled
if do_plot
    figure(1); clf; hold on;
    subplot(1, 2, 1); plot(p, pure_water_density, 'DisplayName', 'Pure Water');
    plot(p, saturated_water_density, 'DisplayName', 'H2 Saturated Water');
    legend; title('Water Density');
    subplot(1, 2, 2); plot(p, pure_gas_density); title('Gas Density');
end

%% Surface densities for calculation
rhoGS = 0.0851; % Density of H2 at surface (kg/m^3)
rhoOS = 9.98207150467e+02; % Density of H2O at surface (kg/m^3)

%% Calculate saturation properties Rs and Rv
if dissolve_gas
    R_s = rhoOS .* X_H2 ./ (rhoGS.*(1-X_H2));
    if do_plot
        figure(1); clf;
        plot(p_sol, R_s); title('Saturated R_s');
    end
end

if vaporize_water
    R_v = rhoGS .* Y_h2o ./ (rhoOS.*(1-Y_h2o));

    if do_plot
        figure(1); clf;
        plot(p_sol, R_v); title('Saturated R_v')
    end
end

%% Compute shrinkage factors for liquid and gas
bO = shrinkage_factor(water_density(T, pure_water_density, X_H2), Y_h2o, rhoOS, rhoGS);
bG = shrinkage_factor(pure_gas_density, Y_h2o, rhoGS, rhoOS);

if do_plot
    figure(1); clf; hold on
    subplot(1, 2, 1)
    plot(p, bG);
    title('Saturated b_g');
    subplot(1, 2, 2)
    plot(p, bO);
    title('Saturated b_o');
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


function writeIMMISC(p, b, mu, title, opts)
% Function to write immiscible properties to a text file
% Inputs:
%   p     - Pressure values (vector)
%   b     - Formation volume factor (vector)
%   mu    - Viscosity values (vector)
%   title - Title/name for the output file (string)
%   opts  - Options structure containing directory path and units

% Construct the output file path
file_path = fullfile(opts.dir, opts.PVDOFile); % Use the filename from opts
fn = fopen(file_path, 'w');

% Check if the file was opened successfully
if fn == -1
    error('Cannot open file: %s', file_path);
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
fprintf('%s written!\n', file_path);
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
file_path = fullfile(opts.dir, opts.PVTOFile); % Use the filename from opts
fn = fopen(file_path, 'w');

nusat = opts.nusat;


% Check if the file was opened successfully
if fn == -1
    error('Cannot open file: %s', file_path);
end

% Write header for the file
fprintf(fn, '-- RS    PRES    FVF      VISC\n');
fprintf(fn, 'PVTO\n');

% Preprocessing for minimum conditions
X_h2_min = 1.0e-8; % Minimum hydrogen saturation
Rs_min = rhoOS * X_h2_min / rhoGS; % Minimum Rs calculation
p_zero = interp1(X_h2_sat, p, X_h2_min, 'pchip', 'extrap');
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
bO_sat = shrinkage_factor(water_density(T, pure_water_density, X_h2_sat), X_h2_sat, rhoOS, rhoGS);

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
        rho_w_dissolved = water_density(T, rho_w_j, X_h2); % Calculate density with dissolved hydrogen
        bO_usat = shrinkage_factor(rho_w_dissolved, X_h2, rhoOS, rhoGS); % Calculate shrinkage factor
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
fprintf('%s written!\n', file_path); % Confirmation message
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
file_path = fullfile(opts.dir, opts.PVTGFile); % Use the filename from opts
fn = fopen(file_path, 'w');

% Check if the file was opened successfully
if fn == -1
    error('Cannot open file: %s', file_path);
end

% Write header for the file
fprintf(fn, '-- PRES       RV       FVF       VISC\n');
fprintf(fn, 'PVTG\n');

% Calculate shrinkage factor once before using it
bG = shrinkage_factor(pure_gas_density, Y_h2o, rhoGS, rhoOS);
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
fprintf('%s written!\n', file_path); % Confirmation message
end



function rho = h2_density_in_water(T)
% Computes the density of hydrogen in water based on temperature.
% Inputs:
% T: Temperature in Celsius (scalar)
% Outputs:
% rho: Density of hydrogen in water (kg/m^3)

% Calculate the partial volume of H2 in water at temperature T
V = h2_partial_volume_in_water(T);
M = 2.016e-3; % Molar mass of H2 (kg/mol)

% Calculate the density of H2
rho = M ./ V; % Density (kg/m^3)
end

function V = h2_partial_volume_in_water(T)
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

function rhoW = water_density(T, pure_water_density, X_h2)
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
vol = X_h2o ./ pure_water_density + X_h2 ./ h2_density_in_water(T);

% Calculate the density of the mixture
rhoW = 1 ./ vol; % Density (kg/m^3)
end

function bO = shrinkage_factor(rhoO, X_h2, rhoOS, rhoGS)
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

function mf = mole_fraction_to_mass_fraction(x, molar_mass_self, molar_mass_other)
% Converts mole fraction to mass fraction.
% Inputs:
% x: Mole fraction (scalar)
% molar_mass_self: Molar mass of the component (kg/mol)
% molar_mass_other: Molar mass of the other component (kg/mol)
% Outputs:
% mf: Mass fraction (scalar)

% Calculate the mass fraction
mf = (x * molar_mass_self) ./ (x .* molar_mass_self + (1 - x) .* molar_mass_other); % Mass fraction
end

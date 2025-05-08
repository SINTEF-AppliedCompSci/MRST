% MRST Example: Setting Up H₂-Brine Fluid Properties for Black-Oil Simulation
% This MRST example demonstrates how to generate black-oil fluid property tables for a 
% hydrogen-brine system under reservoir conditions. The workflow accounts for temperature- 
% and pressure-dependent properties of H₂ and H₂O (optionally including salinity) to compute 
% hydrogen solubility in brine and construct corresponding PVT tables. In addition, the example 
% shows how to generate gas-water flow property (SGOF) tables for black-oil simulation, including 
% gas relative permeability (krG), waqter relative permeability (krO), and capillary pressure at the 
% gas-water contact (pcOG), over a defined saturation range.
%
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
% Parameters for Temperature, Pressure, and Saturation Range
min_temp = 40;                % Minimum temperature in Celsius
max_temp = 60;                % Maximum temperature in Celsius
min_pressure = 1 * atm;       % Minimum pressure in Pa
max_pressure = 45 * barsa();      % Maximum pressure in Pa
nbp = 10;                     % Number of pressure points
nbt = 10;                     % Number of temperature points
ms = 0;                       % Salt molality [mol/kg]
outputDisplay = true;         % Set to true to display generated tables
recompute = true;             % recompuete PVT and SGOF tables

%% Notice on Computational Cost
warning('ComputationalCost:Medium', ...
       ['Please be advised that for large nbp and nbt this example often takes a long time ', ...
        'to run: this script will extract data from https://webbook.nist.gov']);

% Define the target output directory relative to the current directory
outputPathSol = fullfile(mrstOutputDirectory(), 'UHS_PVT', 'H2SolubilityTable');
% Check if the directory exists; if not, create it
if ~exist(outputPathSol, 'dir')
    mkdir(outputPathSol);
end
% Generate H2O Component Table
comp_name = 'H2O';
disp(['Generating component table for: ', comp_name]);
tab_H2O = generateComponentProperties('min_temp',min_temp, 'max_temp',max_temp, 'n_temp', nbt, 'min_press',min_pressure, 'max_press', max_pressure, 'n_press',nbp, 'comp_name', comp_name,'outputDisplay', outputDisplay,'outputPath',outputPathSol);
% Generate H2 Component Table
pause(0.5);  % Ensure smooth execution between commands
comp_name = 'H2';
disp(['Generating component table for: ', comp_name]);
tab_H2 = generateComponentProperties('min_temp',min_temp, 'max_temp',max_temp, 'n_temp', nbt, 'min_press',min_pressure, 'max_press', max_pressure, 'n_press',nbp, 'comp_name', comp_name,'outputDisplay', outputDisplay,'outputPath',outputPathSol);
% Generate Solubility Table
disp('Generating solubility table...');
% [tab_sol, status_sol, file_path_sol] = generateSolubilityTable(min_temp, max_temp, min_pressure, max_pressure, nbp, nbt, ms, outputPath);
tab_sol= generateH2WaterSolubilityTable('min_temp',min_temp, 'max_temp',max_temp, 'n_temp', nbt,'min_press',min_pressure, 'max_press',max_pressure, 'n_press', nbp, 'ms', ms,'outputDisplay', outputDisplay,'outputPath',outputPathSol,'reCompute', recompute);
% Configure and Write Fluid Properties (PVT) Tables
onlyRS = false;
% Define the target output directory relative to the current directory
outputPathPvt = fullfile(mrstOutputDirectory(), 'UHS_PVT', 'PVT_H2STORAGE');
% Check if the directory exists; if not, create it
if ~exist(outputPathPvt, 'dir')
    mkdir(outputPathPvt);
end

if onlyRS
    getFluidH2BrineProps(tab_H2O, tab_H2, tab_sol,'rs', true, 'rv', false,'PVTGFile', 'PVTGH2BRINE', 'PVTOFile', 'PVTOH2BRINE', 'PVTGFile', 'PVTGH2BRINE','outputPath', outputPathPvt, 'reCompute', recompute);
    disp('Writing fluid properties for miscible case with digas and disabled evapoil...');
else
    % generate RSRV
    getFluidH2BrineProps(tab_H2O, tab_H2, tab_sol, 'rs', true, 'rv', true, 'PVTGFile', 'PVTGH2BRINE', 'PVTOFile', 'PVTOH2BRINE','PVDOFile', 'PVDOH2BRINE', 'outputPath', outputPathPvt, 'reCompute', recompute);
    disp('Writing fluid properties for miscible case with digas and evapoil...');
end
% We generates gas-oil flow properties specific to a hydrogen-brine system, for 
% three distinct regions. We exemplify the system with three rock types: caprock, bedrock, and storage rock (aquifer).
% The output includes key properties: gas relative permeability (krG), oil relative permeability (krO),
% and capillary pressure at the gas-oil contact (pcOG). UHS_BENCHMARK_RSRV
%
% Relative permeability is modeled quadratically to capture the flow dynamics between gas and oil phases
% within each rock type, providing a more realistic representation of gas-brine interactions under 
% varying reservoir conditions.
%% Generate SGOF (Gas-Oil Flow) Properties Table
disp('Generating gas-oil flow properties for hydrogen-brine system with three regions');
getFluidH2BrineSGOF('n', 100, 'plot', true, 'outputPath', outputPathPvt, 'fileName', 'SGOF_UHS.txt', ...
                      'units', 'metric', 'nreg', 3, 'sw_imm', [0.1, 0.1, 0.1], ...
                      'sg_imm', [0.05, 0.1, 0.1], 'c_a1', 3.5, 'c_a2', 3.5, 'c_a3', 1.5, ...
                      'Pe', [0.4, 5, 10], 'P_c_max', [9.5e4, 9.5e4, 9.5e4], 'reCompute', recompute);

% Display completion message
disp('H₂-Brine fluid properties setup for black-oil simulation is complete.');

% Force recompute tables
if ~recompute
    disp('You have changed P, T, or ms! Make sure "recompute" is set to true to recalculate tables.');
end

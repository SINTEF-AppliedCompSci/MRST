function tab_comp = generateComponentProperties(varargin)
% MAKE_COMPONENT_TABLE_H2 Generate tables for H2O and H2 fluid properties
% with dynamic input parameters.
%
% SYNOPSIS:
%   tab_comp = MAKE_COMPONENT_TABLE_H2('min_temp', 0, 'max_temp', 100,
%                                      'n_temp', 11, 'min_press', 1e5,
%                                      'max_press', 1e7, 'n_press', 10,
%                                      'comp_name', 'H2O',
%                                      'outputDisplay', false,
%                                      'outputPath', dir)
%
% DESCRIPTION:
%   This function generates a table for the fluid properties of a given
%   component (e.g., H2 or H2O) across a specified temperature and pressure
%   range, using dynamic input parameters for flexibility.
%
% REQUIRED PARAMETERS:
%   'min_temp'    - Minimum temperature in Celsius (default: 0°C)
%   'max_temp'    - Maximum temperature in Celsius (default: 100°C)
%   'n_temp'      - Number of temperature points (default: 11)
%   'min_press'   - Minimum pressure in Pascals (default: 1e5 Pa)
%   'max_press'   - Maximum pressure in Pascals (default: 1e7 Pa)
%   'n_press'     - Number of pressure points (default: 10)
%   'comp_name'   - Name of the component (default: 'H2O')
%
% OPTIONAL PARAMETERS:
%   'outputDisplay'  - If true, displays the generated table. Default: false
%   'outputPath'     - Directory where the output file is saved. Default: current directory.
%   'fileName'       - File name for the generated table.
%
% EXAMPLE USAGE:
%   tab_comp = MAKE_COMPONENT_TABLE_H2('min_temp', 0, 'max_temp', 100, ...
%                                       'n_temp', 11, 'min_press', 1e5, ...
%                                       'max_press', 1e7, 'n_press', 10, ...
%                                       'comp_name', 'H2O', ...
%                                       'outputPath', '/path/to/output', ...
%                                       'fileName', 'water_properties.csv');
%
% COPYRIGHT:
%   Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.
%   This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
%
%   MRST is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   MRST is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MRST. If not, see <http://www.gnu.org/licenses/>.

% Input parameters
opts = struct(...
    'min_temp', 0, ...          % Default: 0°C
    'max_temp', 100, ...        % Default: 100°C
    'n_temp', 11, ...           % Default: 11 points
    'min_press', 1e5, ...       % Default: 1 bar (1e5 Pa)
    'max_press', 1e7, ...       % Default: 100 bar (1e7 Pa)
    'n_press', 10, ...          % Default: 10 points
    'comp_name', 'H2O', ...     % Default: 'H2O'
    'outputDisplay', false, ... % Default: false (auto-generated)
    'outputPath', '', ...       % Default: mrst output directory
    'fileName', '', ...         % Default: see below
    'reCompute', true ...       % Default: true
    );


% Merge user inputs with the defaults
opts = merge_options(opts, varargin{:});

% Rest of the function remains the same...
minTemp = opts.min_temp;
maxTemp = opts.max_temp;
nTemp = opts.n_temp;

delta_temperature = (maxTemp - minTemp) / nTemp;

minPress = opts.min_press;
maxPress = opts.max_press;
nPress = opts.n_press;

fprintf('Component: %s\n', opts.comp_name);
fprintf('Temperature range: %.1f to %.1f °C (%d points)\n', minTemp, maxTemp, nTemp);
fprintf('Pressure range: %.1f to %.1f Pa (%d points)\n', minPress, maxPress, nPress);

if nPress < 2
    error('n_press needs to be at least 2');
end
delta_pressure = (maxPress - minPress)./ (nPress-1);

compName = opts.comp_name;

% Set output filename
if isempty(opts.fileName)
    opts.fileName = sprintf('%svalues_%.1f_to_%.1f_bar_%.1f_to_%.1f_C.csv', ...
    lower(compName), minPress, maxPress, minTemp, maxTemp);
end
% Ensure the output directory exists
if isempty(opts.outputPath)
    opts.outputPath = fullfile(mrstOutputDirectory(), 'UHS_PVT', 'H2SolubilityTable');
end

% Full path to the output file
filePath = fullfile(opts.outputPath, opts.fileName);

% Check if the file already exists
if exist(filePath, 'file')&&~opts.reCompute
    fprintf('File %s already exists. Loading data...\n', opts.fileName);
    try
        tab_comp = readtable(filePath, 'VariableNamingRule', 'preserve');
        if opts.outputDisplay
            disp('Component Table read successfully:');
            display(tab_comp);
        end
        return; % Exit the function as the table already exists
    catch ME
        warning('Failed to read existing file: %s. Proceeding with table generation.', ME.message);
    end
end

if ~exist(opts.outputPath, 'dir')
    mkdir(opts.outputPath);
end

% Open the output file
outFile = fopen(filePath, 'w');
if outFile == -1
    error('Could not open the file: %s', filePath);
end

fprintf(outFile, '# This autogenerated file contains thermodynamical properties of %s.\n', compName);
fprintf(outFile, '# The data has been obtained by querying the NIST Chemistry WebBook https://doi.org/10.18434/T4D303.\n#\n');
fprintf(outFile, '# Concerning temperature and pressure ranges, the following parameters have been used:\n');
fprintf(outFile, '# min temperature = %.1f, max temperature = %.1f, #temperature sampling points = %d\n', ...
    minTemp, maxTemp, nTemp);
fprintf(outFile, '# min pressure = %.1f, max pressure = %.1f, #pressure sampling points = %d\n#\n', ...
    minPress, maxPress, nPress);
fprintf(outFile, '# temperature [°C],     pressure [Pa],   density [kg/m3],  viscosity [Pa.s],   enthalpy [J/kg]\n');

% Get data from NIST
for i = 1:nTemp
    temperature = minTemp + (i-1) * delta_temperature;
    query = struct(...
        'Action', 'Data', ...
        'Wide', 'on', ...
        'ID', 'C7732185', ... % H2O
        'Type', 'IsoTherm', ...
        'Digits', '12', ...
        'PLow', num2str(minPress), ...
        'PHigh', num2str(maxPress), ...
        'PInc', num2str(delta_pressure), ...
        'T', num2str(temperature), ...
        'RefState', 'DEF', ...
        'TUnit', 'C', ...
        'PUnit', 'Pa', ...
        'DUnit', 'kg/m3', ...
        'HUnit', 'kJ/kg', ...
        'WUnit', 'm/s', ...
        'VisUnit', 'uPas', ...
        'STUnit', 'N/m');

    if strcmpi(compName, 'H2')
        query.ID = 'C1333740'; % H2
    end

    queryStr = '';
    fields = fieldnames(query);
    for k = 1:length(fields)
        if k > 1
            queryStr = [queryStr '&'];
        end
        queryStr = [queryStr fields{k} '=' query.(fields{k})];
    end

    try
        response = webread(['https://webbook.nist.gov/cgi/fluid.cgi?' queryStr]);
        responseStr = char(response)';  % Convert uint8 to character array
        % Save the response as a temporary text file and read it
        tempFile = 'temp_data.txt';
        fid = fopen(tempFile, 'w');
        fprintf(fid, '%s', responseStr);
        fclose(fid);
        % Read the table
        data = readtable(tempFile, 'Delimiter', '\t', 'ReadVariableNames', true,'VariableNamingRule', 'preserve');
        % Delete the temporary file
        delete(tempFile);

        % Retrive the necessary data
        phase = data.Phase;
        uniquePhases = unique(phase); % Get unique phase types

        if numel(uniquePhases) > 1
            % Detect phase transitions using string comparison
            phaseTransitionForward = [false; ~strcmp(phase(1:end-1), phase(2:end))];
            phaseTransitionBackward = [~strcmp(phase(1:end-1), phase(2:end)); false];
            isNotPhaseTransition = ~(phaseTransitionForward | phaseTransitionBackward);
        else
            % If there's only one phase, no transitions exist
            phaseTransitionForward = false(size(phase));
            phaseTransitionBackward = false(size(phase));
            isNotPhaseTransition = true(size(phase));
        end



        pressure = data.("Pressure (Pa)")(isNotPhaseTransition);
        density = data.("Density (kg/m3)")(isNotPhaseTransition);
        viscosity = data.("Viscosity (uPa*s)")(isNotPhaseTransition) * 1e-6; % Convert to Pa.s
        enthalpy = data.("Enthalpy (kJ/kg)")(isNotPhaseTransition) * 1000;   % Convert to J/kg

        for j = 1:length(pressure)
            fprintf(outFile, ' %.11e, %.11e, %.11e, %.11e, %.11e\n', ...
                temperature, pressure(j), density(j), viscosity(j), enthalpy(j));
        end
    catch ME
        warning('Failed to retrieve data for T=%.1f°C: %s', temperature, ME.message);
    end
end

fclose(outFile);
fprintf('A file %s has been generated in %s.\n', opts.fileName, opts.outputPath);
% Read the CSV file into a table
try

    tab_comp = readtable(filePath,  'VariableNamingRule', 'preserve');
    disp('Table read successfully:');
    if opts.outputDisplay
        display(tab_comp);
    end
catch ME
    disp('Error occurred while reading the table:');
    disp(ME.message);
    tab_comp = []; % Return empty if reading fails
end
end
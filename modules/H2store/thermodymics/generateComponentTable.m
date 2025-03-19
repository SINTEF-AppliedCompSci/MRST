function [data_table, status, file_path] = generateComponentTable(min_temp, max_temp, min_pressure, max_pressure, nbp, nbt, comp_name, output_dir,varargin)
    % generateComponentTable Generates a component table for the specified component.
    %
    % This function creates a thermodynamic table by running an external Python
    % script with specified parameters for temperature, pressure, and sampling points.
    % The generated file is saved in the specified directory, and the table is
    % read into MATLAB as a table output.
    %
    % Inputs:
    %   min_temp       - Minimum temperature in Celsius.
    %   max_temp       - Maximum temperature in Celsius.
    %   min_pressure   - Minimum pressure in Pascals.
    %   max_pressure   - Maximum pressure in Pascals.
    %   nbp            - Number of pressure sampling points.
    %   nbt            - Number of temperature sampling points.
    %   comp_name      - Name of the component (e.g., 'H2O').
    %   output_dir     - Directory where the generated file will be saved.
    %
    % Outputs:
    %   data_table     - MATLAB table containing the generated data.
    %   status         - Status of the system command execution (0 for success).
    %   file_path      - Full path of the generated CSV file.
    % SEE ALSO:
    %   getFluidH2BrineProps, generateSolubilityTable.

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
% Default options
    opts.outputDisplay = true;  % Default to display output
    opts = merge_options(opts,varargin{:});
    % Convert pressures to barsa for compatibility with filename format
    min_pressure_barsa = convertTo(min_pressure, barsa()); % Minimum pressure in bar
    max_pressure_barsa = convertTo(max_pressure, barsa()); % Maximum pressure in bar

    % Construct the subdirectory path

    if ~isfolder(output_dir)
        mkdir(output_dir); % Create the directory if it doesn't exist

    end

    % Construct the filename based on the defined parameters
    file_name = sprintf('%svalues_%.1f_to_%.1f_bar_%.1f_to_%.1f_C.csv', lower(comp_name), ...
        min_pressure_barsa, max_pressure_barsa, ...
        min_temp, max_temp);

    % Full file path
    file_path = fullfile(output_dir, file_name);

    % Define the command to run the Python script, passing parameters to control output
    command = sprintf('python3 /home/elyes/Documents/Projects/MRST/modules/H2store/thermodymics/python_scripts/make_component_table_H2.py -t1 %.1f -t2 %.1f -nt %d -p1 %.1f -p2 %.1f -np %d -c%s -o %s', ...
        min_temp, max_temp, nbt, min_pressure, max_pressure, nbp, comp_name, file_path);

    % Print the command for debugging
    disp(['Executing command: ', command]);

    % Execute the command using the system function
    [status, result] = system(command);

    % Check the status and display the result
    if status == 0
        disp('Command executed successfully:');
        disp(result);
    else
        disp('Error occurred while executing the command:');
        disp(result);
        return; % Exit the function if the command fails
    end

    % Display path information for verification
    fprintf('File will be saved to: %s\n', file_path);

    % Read the CSV file into a table
    try
        data_table = readtable(file_path);
        disp('Table read successfully:');
        if opts.outputDisplay
            disp(data_table);
        end
    catch ME
        disp('Error occurred while reading the table:');
        disp(ME.message);
        data_table = []; % Return empty if reading fails
    end
end

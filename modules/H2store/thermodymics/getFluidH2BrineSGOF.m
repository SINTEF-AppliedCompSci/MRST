function getFluidH2BrineSGOF(varargin)
    % writeFluidH2BrineSGOF - Write the gas-oil flow properties for hydrogen-brine systems.
    %
    % This function computes and writes the gas-oil flow properties, including
    % gas relative permeability (krG), oil relative permeability (krO), and 
    % pressure at gas-oil contact (pcOG), for a specified number of saturation points.
    % The results are saved to a specified text file, allowing for further analysis 
    % and use in simulations.
    %
    % Syntax:
    %   writeFluidH2BrineSGOF('PropertyName', PropertyValue, ...)
    %
    % Optional Input Parameters:
    %   'n'         - Number of saturation points (default: 100)
    %   'plot'      - Boolean to enable/disable plotting (default: true)
    %   'dir'       - Directory for saving the output file (default: '')
    %   'filename'  - Name of the output file (default: 'SGOF_H2STORAGE.txt')
    %   'units'     - Unit system for output values ('metric' or 'imperial', default: 'metric')
    %   'nreg'      - Number of regions (default: 1)
    %   'sw_imm'    - Array of immiscible water saturations for each region (default: [0.1])
    %   'sg_imm'    - Array of immiscible gas saturations for each region (default: [0.1])
    %   'c_a1'      - Coefficient for oil relative permeability (default: 2)
    %   'c_a2'      - Coefficient for gas relative permeability (default: 2)
    %   'c_a3'      - Coefficient for gas-oil contact pressure (default: 1.5)
    %   'Pe'        - Array of pressure entry values for each region (default: [0.4])
    %   'P_c_max'   - Maximum pressure for each region (default: [9.5e4])
    % SEE ALSO:
    %   generateSolubilityTable, getFluidH2BrineProps.

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
    % Parse input arguments and set defaults
    opt = struct('n', 100, ...
                 'plot', true, ...
                 'dir', '', ...
                 'filename', 'SGOF_H2STORAGE.txt', ...
                 'units', 'metric', ...
                 'nreg', 1, ...  % Default changed to 1 region
                 'sw_imm', 0.1, ...  % Default for single region
                 'sg_imm', 0.1, ...  % Default for single region
                 'c_a1', 2, ...
                 'c_a2', 2, ...
                 'c_a3', 1.5, ...
                 'Pe', 0.4, ...  % Default for single region
                 'P_c_max', 9.5e4);  % Default for single region
             
    opt = merge_options(opt, varargin{:});
    u = unitConversionFactors('si', opt.units);

    % Ensure the output directory exists
    if ~isempty(opt.dir) && ~isfolder(opt.dir)
        mkdir(opt.dir);
    end

    % Prepare figures if plotting is enabled
    if opt.plot
        figure(1); clf; % krO and krG plots
        figure(2); clf; % pcOG plot
    end

    % Open file for writing
    fn = fopen(fullfile(opt.dir, opt.filename), 'w');
    fprintf(fn, '-- SG    KRG    KROG   PCOG\n');
    fprintf(fn, 'SGOF\n');

    % Compute properties for each region
    for reg = 1:opt.nreg
        s0 = opt.sg_imm(reg);  % Single value for immiscible gas saturation
        sat = [0, linspace(s0, 1, opt.n-1)]; % Saturation points
        sg = scale_saturation(sat, s0); % Scaled gas saturation
        so = scale_saturation(1 - sat, opt.sw_imm(reg)); % Scaled oil saturation
        krg = sg.^opt.c_a2; % Gas relative permeability
        kro = so.^opt.c_a1; % Oil relative permeability
        pc_og_bar = opt.Pe(reg) * so.^(1 - opt.c_a3); % Pressure at gas-oil contact

        pcmax = opt.P_c_max(reg);
        pc_og = pcmax * erf((pc_og_bar / pcmax) * (sqrt(pi) / 2)); % Adjusted pressure

        % Plotting results if enabled
        if opt.plot
            figure(1);
            subplot(1, 2, 1); hold on;
            title('krO');
            plot(1 - sat, kro);
            subplot(1, 2, 2); hold on;
            title('krG');
            plot(sat, krg);
            figure(2);
            subplot(1, opt.nreg, reg); hold on;
            title('pcOG');
            plot(sat, pc_og_bar, 'o');
            plot(sat, pc_og);
        end
        
        % Write computed properties to file
        for i = 1:opt.n
            fprintf(fn, '%1.8f %1.8f %1.8f %1.8f\n', sat(i), krg(i), kro(i), pc_og(i) * u.press);
        end
        fprintf(fn, '/\n');
    end
    
    % Display completion message
    disp([opt.filename ' written!']);
    fclose(fn); % Close the output file
end

function s_scaled = scale_saturation(s, s_imm)
    % scale_saturation - Scale saturation values.
    %
    % This function scales the saturation values to a normalized range
    % between 0 and 1 based on the immiscible saturation limit.
    %
    % Syntax:
    %   s_scaled = scale_saturation(s, s_imm)
    %
    % Input Parameters:
    %   s (array)   - Saturation values to be scaled.
    %   s_imm (scalar) - Immiscible saturation threshold.
    %
    % Output:
    %   s_scaled (array) - Scaled saturation values, clipped to a minimum of 0.
    
    s_scaled = max((s - s_imm) / (1 - s_imm), 0);
end

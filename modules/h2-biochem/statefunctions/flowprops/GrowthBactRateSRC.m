classdef GrowthBactRateSRC < StateFunction
    % Bacterial growth rate computation for compositional simulations
    %
    % SYNOPSIS:
    %   gr = GrowthBactRateSRC(model, 'property1', value1, ...)
    %
    % DESCRIPTION:
    %   Computes the bacterial growth rate in each grid cell based on:
    %   - H2 and CO2 concentrations (Monod kinetics) for Methacenic archae
    %   - H2 and SO4(2-) concentrations (Monod kinetics) for Sulfate reducing bacteria
    %   - Bacterial population density
    %   - Liquid phase saturation and density
    %   - Pore volume
    %
    % REQUIRED PARAMETERS:
    %   model - Reservoir model with bacterial modeling enabled
    %
    % OPTIONAL PARAMETERS:
    %   None
    %
    % RETURNS:
    %   Class instance ready for use in simulation
    %
    % SEE ALSO:
    %   CompositionalModel, EquationsCompositional

    properties
        % No additional properties needed
    end

    methods
        function gp = GrowthBactRateSRC(model, varargin)
            % Constructor for bacterial growth rate calculator

            % Initialize base class
            gp@StateFunction(model, varargin{:});

            % Declare dependencies
            gp = gp.dependsOn({'x'}, 'state');        % Mole fractions
            gp = gp.dependsOn({'s'}, 'state');        % Saturations
            gp = gp.dependsOn({'nbact'}, 'state');    % Bacterial concentration
            gp = gp.dependsOn({'PoreVolume', 'Density'}, 'PVTPropertyFunctions');

            gp.label = '\Psi_{growth}'; % LaTeX-style label
        end

        function Psigrowth = evaluateOnDomain(prop, model, state)
            % Compute bacterial growth rate in each grid cell
            %
            % PARAMETERS:
            %   prop  - Property function instance
            %   model - Reservoir model instance
            %   state - State struct containing fields
            %
            % RETURNS:
            %   Psigrowth - Bacterial growth rate per cell [1/s]

            % Initialize with zero growth rate
            Psigrowth = 0;

            % Get component names and indices
            rm = model.ReservoirModel;
            namecp = rm.getComponentNames();
            idx_H2 = find(strcmpi(namecp, 'H2'), 1);     % Case-insensitive search
            idx_CO2 = find(strcmpi(namecp, 'CO2'), 1);   % Case-insensitive search

            % Check if bacterial modeling is active and components exist
            if ~(rm.bacteriamodel && rm.liquidPhase && ~isempty(idx_H2) && ~isempty(idx_CO2))
                return;
            end

            % Get required properties
            pv = rm.PVTPropertyFunctions.get(rm, state, 'PoreVolume');
            rho = rm.PVTPropertyFunctions.get(rm, state, 'Density');
            s = rm.getProp(state, 's');
            nbact = rm.getProp(state, 'nbact');
            x = rm.getProp(state, 'x');
            L_ix = rm.getLiquidIndex();

            % Extract liquid phase properties
            if iscell(x)
                xH2 = x{idx_H2};
                xCO2 = x{idx_CO2};
                sL = s{L_ix};
                rhoL = rho{L_ix};
            else
                xH2 = x(:, idx_H2);
                xCO2 = x(:, idx_CO2);
                sL = s(:, L_ix);
                rhoL = rho(:, L_ix);
            end

            % Calculate effective volume with safeguards
            if iscell(sL)
                Voln = max(sL{1}, 1.0e-8) .* rhoL{1};
            else
                Voln = max(sL, 1.0e-8) .* rhoL;
            end
            Voln = max(Voln, 1.0e-8);

            % Get growth parameters
            alphaH2 = rm.alphaH2;
            alphaCO2 = rm.alphaCO2;
            Psigrowthmax = rm.Psigrowthmax;

            % Calculate Monod terms for H2 and CO2
            axH2 = xH2 ./ (alphaH2 + xH2);
            axCO2 = xCO2 ./ (alphaCO2 + xCO2);

            % Compute growth rate
            Psigrowth = pv .* Psigrowthmax .* axH2 .* axCO2 .* nbact .* Voln;
        end
    end
end

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
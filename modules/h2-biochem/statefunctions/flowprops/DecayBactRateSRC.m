classdef DecayBactRateSRC < StateFunction
    % Bacterial decay rate computation for compositional simulations
    %
    % SYNOPSIS:
    %   decay = DecayBactRateSRC(model, 'property1', value1, ...)
    %
    % DESCRIPTION:
    %   Computes the bacterial decay rate in each grid cell, accounting for:
    %   - Bacterial concentration (nbact)
    %   - Liquid phase saturation
    %   - Pore volume
    %   - Phase densities
    %   - Presence of required components (H2 and CO2)
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
        function gp = DecayBactRateSRC(model, varargin)
            % Constructor for bacterial decay rate calculator
            gp@StateFunction(model, varargin{:});

            % Define dependencies
            gp = gp.dependsOn({'nbact'}, 'state');          % Bacterial concentration
            gp = gp.dependsOn({'s'}, 'state');             % Phase saturations
            gp = gp.dependsOn({'PoreVolume', 'Density'}, 'PVTPropertyFunctions');

            % Set label for output
            gp.label = 'Psi_{decay}';
        end

        function Psidecay = evaluateOnDomain(prop, model, state)
            % Compute bacterial decay rate in each grid cell
            %
            % PARAMETERS:
            %   prop  - Property function instance
            %   model - Reservoir model instance
            %   state - State struct with fields
            %
            % RETURNS:
            %   Psidecay - Bacterial decay rate per cell [1/s]

            % Initialize with zeros
            Psidecay = 0;

            % Get model parameters
            rm = model.ReservoirModel;
            bbact = rm.b_bact;
            nbMax = rm.nbactMax;

            % Check if bacterial modeling is active
            if ~(rm.bacteriamodel && rm.liquidPhase)
                return;
            end

            % Get component names and indices
            namecp = rm.getComponentNames();
            idx_H2 = find(strcmpi(namecp, 'H2'), 1);
            idx_CO2 = find(strcmpi(namecp, 'CO2'), 1);

            % Validate required components
            if isempty(idx_H2) || isempty(idx_CO2)
                return;
            end

            % Get required state variables
            pv = rm.PVTPropertyFunctions.get(rm, state, 'PoreVolume');
            rho = rm.PVTPropertyFunctions.get(rm, state, 'Density');
            s = rm.getProp(state, 's');
            nbact = rm.getProp(state, 'nbact');
            L_ix = rm.getLiquidIndex();

            % Extract liquid phase properties
            if iscell(s)
                sL = s{L_ix};
                rhoL = rho{L_ix};
            else
                sL = s(:, L_ix);
                rhoL = rho(:, L_ix);
            end

            % Calculate effective volume with safeguards
            if iscell(sL)
                Voln = max(sL{1}, 1.0e-8) .* rhoL{1};
            else
                Voln = max(sL, 1.0e-8) .* rhoL;
            end

            % Compute decay rate
            Psidecay = pv .* bbact .* nbact .* (nbact .* Voln);

            % Handle negative bacterial concentrations
            Psidecay(nbact < 0) = -Psidecay(nbact < 0);
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
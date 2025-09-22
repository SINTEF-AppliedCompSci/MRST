classdef DynamicFlowPoreVolume < PoreVolume
    % Effective pore-volume accounting for pressure and bacterial effects
    %
    % SYNOPSIS:
    %   pv = DynamicFlowPoreVolume(model, 'property1', value1, ...)
    %
    % DESCRIPTION:
    %   Computes effective pore volume accounting for:
    %   - Rock compressibility effects (via pvMultR)
    %   - Bacterial concentration effects (if bacteriamodel enabled)
    %
    % REQUIRED PARAMETERS:
    %   model - Reservoir model with fluid.pvMultR defined
    %
    % OPTIONAL PARAMETERS:
    %   None
    %
    % RETURNS:
    %   Class instance for dynamic pore volume calculation
    %
    % SEE ALSO:
    %   PoreVolume, ReservoirModel

    properties
        % No additional properties needed
    end

    methods
        function gp = DynamicFlowPoreVolume(model, varargin)
            % Constructor for dynamic pore volume calculator

            % Initialize base class
            gp@PoreVolume(model, varargin{:});

            % Declare dependencies based on model configuration
            if model.bacteriamodel
                gp = gp.dependsOn('nbact', 'state');
            end
            gp = gp.dependsOn('pressure', 'state');

            % Validate required fluid property exists
            assert(isfield(model.fluid, 'pvMultR'), ...
                'Rock compressibility function (pvMultR) missing from fluid.');
        end

        function pv = evaluateOnDomain(prop, model, state)
            % Compute effective pore volume with pressure/bacteria effects
            %
            % PARAMETERS:
            %   prop  - Property function instance
            %   model - Reservoir model instance
            %   state - State struct containing fields
            %
            % RETURNS:
            %   pv - Effective pore volume values

            % Get base pore volume from parent class
            pv = evaluateOnDomain@PoreVolume(prop, model, state);

            % Get current pressure
            p = model.getProp(state, 'pressure');

            % Apply rock compressibility multiplier
            if model.bacteriamodel
                % Include bacterial concentration effect if enabled
                nbact = model.getProp(state, 'nbact');
                pvMult = prop.evaluateFluid(model, 'pvMultR', p, nbact);
            else
                % Pressure-only effect
                pvMult = prop.evaluateFluid(model, 'pvMultR', p);
            end

            % Apply multiplier
            pv = pv .* pvMult;
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
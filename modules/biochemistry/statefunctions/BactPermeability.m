classdef BactPermeability < StateFunction
    % Pressure- and bacteria-dependent permeability modifier
    %
    % SYNOPSIS:
    %   perm = BactPermeability(model)
    %
    % DESCRIPTION:
    %   Computes permeability modification based on:
    %   - Pressure effects (if dynamicFlowTrans enabled)
    %   - Bacterial concentration effects (if bacteriamodel enabled)
    %
    % REQUIRED PARAMETERS:
    %   model - Reservoir model with appropriate settings
    %
    % OPTIONAL PARAMETERS:
    %   None (configure through model settings)
    %
    % RETURNS:
    %   Class instance for permeability calculation
    %
    % SEE ALSO:
    %   ReservoirModel, Rock

    properties
        % No additional properties needed
    end

    methods
        function perm = BactPermeability(model)
            % Constructor for bacteria-permeability relationship
            perm@StateFunction(model);

            % Declare dependencies based on model configuration
            if model.dynamicFlowTrans
                perm = perm.dependsOn('pressure', 'state');
                if model.bacteriamodel
                    perm = perm.dependsOn('nbact', 'state');
                end
            end
            perm.label = 'K';  % Permeability label
        end

        function perm = evaluateOnDomain(prop, model, state)
            % Evaluate permeability modification
            %
            % PARAMETERS:
            %   prop  - Property function instance
            %   model - Reservoir model instance
            %   state - State struct containing fields
            %
            % RETURNS:
            %   perm - Modified permeability values

            % Start with base rock permeability
            perm = model.rock.perm;

            % Apply modifications if enabled
            if model.dynamicFlowTrans
                if model.bacteriamodel
                    % Get both pressure and bacteria concentration
                    [p, nbact] = model.getProps(state, 'pressure', 'nbact');
                    perm = perm(p, nbact);  % Apply both modifications
                else
                    % Pressure-only modification
                    p = model.getProps(state, 'pressure');
                    perm = perm(p, 0);  % Zero bacterial effect
                end
            end
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
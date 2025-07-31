classdef BactPorosity < StateFunction
    % Pressure- and bacteria-dependent porosity modifier
    %
    % SYNOPSIS:
    %   poro = BactPorosity(model)
    %
    % DESCRIPTION:
    %   Computes porosity modification based on:
    %   - Pressure effects (if dynamicFlowPv enabled)
    %   - Bacterial concentration effects (if bacteriamodel enabled)
    %
    % REQUIRED PARAMETERS:
    %   model - Reservoir model with appropriate settings
    %
    % OPTIONAL PARAMETERS:
    %   None (configure through model settings)
    %
    % RETURNS:
    %   Class instance for porosity calculation
    %
    % SEE ALSO:
    %   ReservoirModel, Rock

    properties
        % No additional properties needed
    end
   
    methods       
        function poro = BactPorosity(model)
            % Constructor for bacteria-porosity relationship
            poro@StateFunction(model);
            
            % Declare dependencies based on model configuration
            if model.dynamicFlowPv
                if model.bacteriamodel
                    poro = poro.dependsOn({'pressure', 'nbact'}, 'state');
                else
                    poro = poro.dependsOn('pressure', 'state');
                end
            end
            poro.label = 'Phi';
        end
       
        function poro = evaluateOnDomain(prop, model, state)
            % Evaluate porosity modification
            %
            % PARAMETERS:
            %   prop  - Property function instance
            %   model - Reservoir model instance
            %   state - State struct containing fields
            %
            % RETURNS:
            %   poro - Modified porosity values
            
            % Start with base rock porosity
            poro = model.rock.poro;
            
            % Apply modifications if enabled
            if model.dynamicFlowPv
                if model.bacteriamodel
                    % Get both pressure and bacteria concentration
                    [p, nbact] = model.getProps(state, 'pressure', 'nbact');
                    poro = poro(p, nbact); % Apply both modifications
                else
                    % Pressure-only modification
                    p = model.getProps(state, 'pressure');
                    poro = poro(p, 0); % Zero bacterial effect
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
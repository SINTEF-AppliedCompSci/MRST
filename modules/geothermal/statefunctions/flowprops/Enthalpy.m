classdef Enthalpy < StateFunction
% State function for phase enthalpy in geothermal models
%
%   Computes the specific enthalpy for each phase in the reservoir model.
%   In temperature-based thermal formulations, enthalpy is calculated from
%   internal energy and pressure/density. Used in energy balance equations
%   and for post-processing thermal results.
%
%   Usage:
%      h = Enthalpy(model)
%   where 'model' is a ReservoirModel with a thermal formulation.
%
%   In 'temperature' formulation, enthalpy is computed as:
%      h = u + p/rho
%   where u is the phase internal energy, p is phase pressure, and rho is phase density.
%
%   Dependencies:
%      - PhasePressures
%      - Density
%      - PhaseInternalEnergy (from FlowPropertyFunctions)
%
%   See also: PhaseInternalEnergy, ReservoirModel

    properties
        % No additional properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = Enthalpy(model, varargin)
        % Constructor for Enthalpy state function
        %   Sets up dependencies based on the thermal formulation.
        %
        %   Inputs:
        %      model    - ReservoirModel object
        %      varargin - Additional arguments passed to StateFunction
        %
        %   The dependencies are set according to the model's thermal formulation.
            gp@StateFunction(model, varargin{:});
            switch model.thermalFormulation
                case 'enthalpy'
                    % Enthalpy is the primary variable; no dependencies needed
                case 'temperature'
                    % Temperature is the primary variable; enthalpy is computed
                    gp = gp.dependsOn({'PhasePressures', 'Density'});
                    gp = gp.dependsOn('PhaseInternalEnergy', 'FlowPropertyFunctions');
            end
            gp.label = 'h_\alpha';
        end
        
        %-----------------------------------------------------------------%
        function h = evaluateOnDomain(prop, model, state)
        %evaluateOnDomain   Compute phase enthalpy for each phase
        %   h = evaluateOnDomain(prop, model, state) returns a cell array of
        %   enthalpy values for each phase, based on the current state and
        %   thermal formulation.
        %
        %   Inputs:
        %      prop  - This Enthalpy state function
        %      model - ReservoirModel object
        %      state - Current simulation state
        %
        %   Output:
        %      h     - Cell array of enthalpy values for each phase
            switch model.thermalFormulation
                case 'enthalpy'
                    error('This line should not be reached - something is wrong');
                case 'temperature'
                    [p, rho] = prop.getEvaluatedDependencies(state, 'PhasePressures', 'Density');
                    u = model.getProps(state, 'PhaseInternalEnergy');
                    nph = model.getNumberOfPhases();
                    h = cell(1, nph);
                    for i = 1:nph
                        % Enthalpy: h = u + p/rho for each phase
                        h{i} = u{i} + p{i}./rho{i};
                    end
            end
        end
    end

end

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
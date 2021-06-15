classdef OilWaterSurfactantModel < ThreePhaseSurfactantPolymerModel
    %
    %
    % SYNOPSIS:
    %   model = OilWaterSurfactantModel(G, rock, fluid, varargin)
    %
    % DESCRIPTION: 
    %   Fully implicit model for an oil water system with surfactant. All
    %   the equations are solved implicitly. A description of the surfactant model
    %   that is implemented can be found in the directory ad-eor/docs .
    %
    % PARAMETERS:
    %   G        - Grid
    %   rock     - Rock structure
    %   fluid    - Fluid structure
    %   varargin - optional parameter
    %
    % RETURNS:
    %   class instance
    %
    % EXAMPLE:
    %
    % SEE ALSO: equationsOilWaterSurfactant,
    %

    properties
    end

    methods

        function model = OilWaterSurfactantModel(G, rock, fluid, varargin)

            model = model@ThreePhaseSurfactantPolymerModel(G, rock, fluid);
            model.surfactant = true;
            model.gas = false;
            model.water = true;
            model.oil = true;
            model.disgas = false;
            model.vapoil = false;
            model.polymer = false;
            
            model = merge_options(model, varargin{:});
            model.operators.veloc   = computeVelocTPFA(G, model.operators.internalConn);
            model.operators.sqVeloc = computeSqVelocTPFA(G, model.operators.internalConn);
        end

        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsOilWaterSurfactant(state0, state, ...
                                                           model, dt, ...
                                                           drivingForces, ...
                                                           varargin{:});
        end
    end
end


%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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


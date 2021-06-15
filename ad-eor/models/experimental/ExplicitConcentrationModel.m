classdef ExplicitConcentrationModel < OilWaterSurfactantBaseModel
%
%
% SYNOPSIS:
%   model = ExplicitConcentrationModel(G, rock, fluid, varargin)
%
% DESCRIPTION: Model used to solve for the surfactant concentration using the
% transport equation for an oil-water system with surfactant. The time
% discretization is explicit.
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
% SEE ALSO: ImplicitExplicitOilWaterSurfactantModel, equationsSurfactantTransport
%


    methods
        function model = ExplicitConcentrationModel(G, rock, fluid, varargin)
            model = model@OilWaterSurfactantBaseModel(G, rock, fluid, varargin{:});
            model.wellVarNames = {'qWSft', 'bhp'};
        end


        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsSurfactantTransport(state0, state, model, dt, drivingForces, varargin{:});
        end

        function [state, report] = updateState(model, state, problem, dx, drivingForces)

            [state, report] = updateState@OilWaterSurfactantBaseModel(model, state, problem,  dx, ...
                                                              drivingForces);
            % cap the concentration (only if implicit solver for concentration)
            c = model.getProp(state, 'surfactant');
            state = model.setProp(state, 'surfactant', max(c, 0) );
            
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@OilWaterSurfactantBaseModel(model, state0, state, ...
                                                              dt, drivingForces);
            state = updateAdsorption(state0, state, model);
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

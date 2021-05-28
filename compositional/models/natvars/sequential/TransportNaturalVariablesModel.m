classdef TransportNaturalVariablesModel < NaturalVariablesCompositionalModel
    properties
        staticUpwind
        upwindType
    end
    
    methods
        function model = TransportNaturalVariablesModel(G, rock, fluid, compFluid, varargin)
            
            model = model@NaturalVariablesCompositionalModel(G, rock, fluid, compFluid);
            model.staticUpwind = false;
            model.upwindType  = 'potential';
            model = merge_options(model, varargin{:});
            model.allowLargeSaturations = true;
            model.maxPhaseChangesNonLinear = inf;
        end

        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = transportEquationsNaturalVars(state0, state, model, dt, ...
                            drivingForces, varargin{:});
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

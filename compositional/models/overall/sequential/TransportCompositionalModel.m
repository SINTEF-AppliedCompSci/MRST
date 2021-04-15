classdef TransportCompositionalModel < ThreePhaseCompositionalModel
    % Two phase oil/water system without dissolution
    properties
        conserveWater
        staticUpwind
        upwindType
    end
    
    methods
        function model = TransportCompositionalModel(G, rock, fluid, compfluid, varargin)
            
            model = model@ThreePhaseCompositionalModel(G, rock, fluid, compfluid);
            
            model.conserveWater = true;
            model.staticUpwind = false;
            model.upwindType  = 'potential';

            model = merge_options(model, varargin{:});
            if model.water
                model.saturationVarNames = {'sw', 'so', 'sg'};
            else
                model.saturationVarNames = {'so', 'sg'};
            end
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = transportEquationCompositional(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            'solveForWater', model.conserveWater, ...
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

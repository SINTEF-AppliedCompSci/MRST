classdef TransportOilWaterModel < TwoPhaseOilWaterModel
    % Two phase oil/water system without dissolution
    properties
        conserveWater
        conserveOil
        upwindType
    end
    
    methods
        function model = TransportOilWaterModel(G, rock, fluid, varargin)
            
            model = model@TwoPhaseOilWaterModel(G, rock, fluid);
            
            model.conserveWater = true;
            model.conserveOil   = false;
            
            model.upwindType = 'potential';
            model.useCNVConvergence = false;
            model.nonlinearTolerance = 1e-3;
            
            model = merge_options(model, varargin{:});
            
            assert(~(model.conserveWater && model.conserveOil), ... 
                            'Sequential form only conserves n-1 phases');
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = transportEquationOilWater(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            'solveForOil',   model.conserveOil, ...
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

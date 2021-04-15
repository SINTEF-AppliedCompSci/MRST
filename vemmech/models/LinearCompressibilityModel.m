classdef LinearCompressibilityModel < ReservoirModel
    properties
        thermal = false;
    end
    
    methods
        function model = LinearCompressibilityModel(G, rock, fluid, varargin)
            model = model@ReservoirModel(G, rock, fluid);
            % Only enable water
            model.oil = false;
            model.gas = false;
            model.water = true;
            model.useCNVConvergence = false;
            
            model = merge_options(model, varargin{:});
            
            model.stepFunctionIsLinear = true;
        end

        
        function [state, report] = stepFunction(model, state, state0, dt, drivingForces, linsolve, nonlinsolve, iteration, varargin)

        T = model.operators.T_all;
        
        solfn = @(A, b) linsolve.solveLinearSystem(A, b);
        state = lincompTPFA(dt, state0, model.G, T, model.operators.pv, model.fluid, model.rock,...
            'src', drivingForces.src,...
            'bc', drivingForces.bc, ...
            'Wells', drivingForces.W, ...
            'LinSolve', solfn,...
            'use_trans', true);
        
        if model.thermal
            assert(all(isfinite(model.rock.lambdaR)), ...
                ['Thermal conductivity is not finite and thermal is enabled.'...
                ' Aborting.']);
            rock_heat = struct('perm', model.rock.lambdaR);
            T_r = computeTrans(model.G, rock_heat);
            state = linearTransport(dt, state, model.G, T_r, model.operators.pv, model.fluid, model.rock,...
                            'src', drivingForces.src,...
                            'bc', drivingForces.bc, ...
                            'Wells', drivingForces.W ...
                            );
        end
            
        report = model.makeStepReport(...
                        'Failure',      false, ...
                        'Converged',    true, ...
                        'Residuals',    0, ...
                        'ResidualsConverged', true);
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

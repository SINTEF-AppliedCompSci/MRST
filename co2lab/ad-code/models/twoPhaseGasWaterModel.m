classdef twoPhaseGasWaterModel < ReservoirModel

    properties
        % (constant) temperature field
        t
        name
    end
    
    % ============================================================================
    methods
        % ------------------------------------------------------------------------
        function model = twoPhaseGasWaterModel(G, rock, fluid, tsurf, tgrad, varargin)
           
            model = model@ReservoirModel(G); 
           
            %model.G     = G;
            model.fluid = fluid;
            model.oil   = false;
            model.gas   = true;
            model.water = true;
            model.t     = computeTemperatureField(G, tsurf, tgrad);
            model.name  = 'GasWater_2ph';
            model.saturationVarNames = {'sw', 'sg'}; % @@ Design: ideally, we
                                                     % should not have to set
                                                     % both this _and_ the
                                                     % 'water', 'gas' and 'oil'
                                               % flags above. Check with
                                               % maintainer of parent class.
            model.wellVarNames = {'qWs', 'qGs', 'bhp'};
            model.gravity = gravity;
            model       = model.setupOperators(G, rock);
        end

    end
    
    methods %(Access = protected)

        % ------------------------------------------------------------------------
        
        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, ...
                                                        varargin)
            [problem, state] = equationsGasWater(model, state0, state , dt , ...
                                                 drivingForces             , ...
                                                 varargin{:});
        end

        % ------------------------------------------------------------------------
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces) %#ok
            
           [state, report] = updateState@ReservoirModel(model, state, problem, dx, ...
                                                        drivingForces);
           state.sGmax = max(state.sGmax,state.s(:,2));
           state.sGmax = min(1,state.sGmax);
           state.sGmax = max(0,state.sGmax);
        end
    end
    
end
% ============================================================================

function res = truncate(m, val, M)
% ensure all components in 'val' are in the range [m, M]
    if ~isempty(val)
        res = max(m, min(val, M)); 
    else
        res = [];
    end
end

function t = computeTemperatureField(G, tsurf, tgrad)
    t = tsurf * ones(G.cells.num, 1) + G.cells.centroids(:,3) * tgrad / 1000;
end

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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


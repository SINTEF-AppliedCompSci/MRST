classdef TwoPhaseWaterGasModel < ThreePhaseBlackOilModel
    % Two-phase gas and water model
    properties
        % (constant) temperature field
        t
        name
    end
    
    % ============================================================================
    methods
        % ------------------------------------------------------------------------
        function model = TwoPhaseWaterGasModel(G, rock, fluid, tsurf, tgrad, varargin)
           
            model = model@ThreePhaseBlackOilModel(G, rock, fluid); 

            model.oil   = false;
            model.gas   = true;
            model.water = true;
            if nargin < 4
                [tsurf, tgrad] = deal(nan);
            end
            model.t     = computeTemperatureField(G, tsurf, tgrad);
            model.name  = 'GasWater_2ph';
            model.gravity = gravity;
            model = merge_options(model, varargin{:});
        end

    end
    
    methods

        % ------------------------------------------------------------------------
        
        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, ...
                                                        varargin)
            [problem, state] = equationsWaterGas(model, state0, state , dt , ...
                                                 drivingForces             , ...
                                                 varargin{:});
        end

        % ------------------------------------------------------------------------
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            
           [state, report] = updateState@ThreePhaseBlackOilModel(model, state, problem, dx, ...
                                                        drivingForces);
           sG = model.getProp(state, 'sG');
           if isfield(state, 'sGmax')
               state.sGmax = max(state.sGmax, sG);
           else
               state.sGmax = sG;
           end
           state.sGmax = min(1,state.sGmax);
           state.sGmax = max(0,state.sGmax);
        end
    end
    
end
% ============================================================================

function t = computeTemperatureField(G, tsurf, tgrad)
    t = tsurf * ones(G.cells.num, 1) + G.cells.centroids(:,3) * tgrad / 1000;
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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


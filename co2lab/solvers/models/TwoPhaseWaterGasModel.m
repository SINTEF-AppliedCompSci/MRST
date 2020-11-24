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
        
        % --------------------------------------------------------------------%
        function scaling = getScalingFactorsCPR(model, problem, names, solver)
            % Get approximate, impes-like pressure scaling factors
            nNames = numel(names);

            scaling = cell(nNames, 1);
            handled = false(nNames, 1);

            % Take averaged pressure for scaling factors
            state = problem.state;
            fluid = model.fluid;
            if (isprop(solver, 'trueIMPES') || isfield(solver, 'trueIMPES')) && solver.trueIMPES
                % Rigorous pressure equation (requires lots of evaluations)
                p = state.pressure;
            else
                % Very simple scaling factors, uniform over grid
                p = mean(state.pressure);
            end
            for iter = 1:nNames
                nm = lower(names{iter});
                switch nm
                    case 'water'
                        bW = fluid.bW(p);
                        s = 1./bW;
                    case 'gas'
                        bG = fluid.bG(p);
                        s = 1./bG;
                    otherwise
                        continue
                end
                sub = strcmpi(problem.equationNames, nm);

                scaling{iter} = s;
                handled(sub) = true;
            end
            if ~all(handled)
                % Get rest of scaling factors from parent class
                other = getScalingFactorsCPR@ReservoirModel(model, problem, names(~handled));
                [scaling{~handled}] = other{:};
            end
        end

    end
    
end
% ============================================================================

function t = computeTemperatureField(G, tsurf, tgrad)
    t = tsurf * ones(G.cells.num, 1) + G.cells.centroids(:,3) * tgrad / 1000;
end

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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


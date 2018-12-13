classdef TwoPhaseOilWaterModel < ThreePhaseBlackOilModel
    % Two phase oil/water system without dissolution
    properties

    end

    methods
        function model = TwoPhaseOilWaterModel(G, rock, fluid, varargin)
            model = model@ThreePhaseBlackOilModel(G, rock, fluid);

            % This is the model parameters for oil/water
            model.oil = true;
            model.gas = false;
            model.water = true;

            model = merge_options(model, varargin{:});
        end

        % --------------------------------------------------------------------%
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            if 0
            [problem, state] = equationsOilWater(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});
            else
                opt = struct('Verbose',     mrstVerbose,...
                            'reverseMode', false,...
                            'resOnly',     false,...
                            'iteration',   -1);
                opt = merge_options(opt, varargin{:});


                % Define primary variables
                if opt.resOnly
                    [dyn_state, primaryVars] = model.getDynamicState(state);
                    [dyn_state0, ~]          = model.getDynamicState(state0);
                else
                    if ~opt.reverseMode
                        [dyn_state, primaryVars] = model.getForwardDynamicState(state);
                        [dyn_state0, ~]          = model.getDynamicState(state0);
                    else
                        [dyn_state, ~]            = model.getDynamicState(state);
                        [dyn_state0, primaryVars] = model.getReverseDynamicState(state0);
                    end
                end

                [eqs, names, types] = equationsOilWaterDynamicState(dyn_state0, ...
                                                                  dyn_state, ...
                                                                  model, dt, ...
                                                                  drivingForces);

                dissolved = model.getDissolutionMatrix(dyn_state.rs, dyn_state.rv);
                % Add in and setup well equations

                ws_dyn   = dyn_state.wellSol;
                wellVars = ws_dyn.dynamicVariables(1:end-1);
                wellMap  = ws_dyn.wellmap;

                p   = dyn_state.pressure;
                mob = dyn_state.FlowProps.Mobility;
                rho = dyn_state.FlowProps.Density;
                
                [eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, ...
                                                                  names, types, ...
                                                                  state0.wellSol, ...
                                                                  state.wellSol, ...
                                                                  wellVars, ...
                                                                  wellMap, p, ...
                                                                  mob, rho, ...
                                                                  dissolved, ...
                                                                  {}, dt, opt);

                state.FlowProps = dyn_state.FlowProps.reduce();

                problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

            end
        end
        
        function [dyn_state, primaryVariables] = getForwardDynamicState(model, state)
            [dyn_state, primaryVariables] = setupDynamicStateOilWater(model, ...
                                                              state, true, false);
        end

        function [dyn_state, primaryVariables] = getDynamicState(model, state)
            [dyn_state, primaryVariables] = setupDynamicStateOilWater(model, ...
                                                              state, false);
        end

        function [dyn_state, primaryVariables] = getReverseDynamicState(model, state)
            [dyn_state, primaryVariables] = setupDynamicStateOilWater(model, ...
                                                              state, true, true);
        end
    end
end

%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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
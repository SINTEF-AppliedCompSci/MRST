classdef MechFluidFixedStressSplitModel < MechFluidSplitModel
%
%
% SYNOPSIS:
%   MechFluidFixedStressSplitModel(G, rock, fluid, mech_problem, varargin)
%
% DESCRIPTION:
%   Model for the fixed stress splitting method. The model contains a
%   mechanic and fluid model, see MechFluidSplitModel. Different fluid
%   models can be used, see setupFluidModel member function.
%
%   This model essentially overwrites the member function stepFunction.
%   There, the mechanic equations are solved for a given pore (fluid)
%   pressure. Then, a strain-pressure relation is established (see
%   uvcomputeMechTerm member function) assuming that the total stress is
%   constant. Finally, the fluid equation is solved.
%
% PARAMETERS:
%   G            - grid structure
%   rock         - rock structure
%   fluid        - fluid structure
%   mech_problem - structure containing mechanical parameters
%   varargin     - 
%
% RETURNS:
%   class instance
%
% EXAMPLE:
%   run2DCase, runNorneExample
%
% SEE ALSO:
%   MechFluidSplitModel, BlackOilFixedStressFluidModel,
%   OilWaterFixedStressFluidModel, WaterFixedStressFluidModel

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

    methods
        function model = MechFluidFixedStressSplitModel(G, rock, fluid, ...
                                                        mech_problem, varargin)

            model = model@MechFluidSplitModel(G, rock, fluid, mech_problem, ...
                                         varargin{:});
        end

        function fluidModel = setupFluidModel(model, rock, fluid, fluidModelType, ...
                                                     varargin)
            switch fluidModelType
              case 'water'
                fluidModel = WaterFixedStressFluidModel(model.G, rock, fluid, varargin{:});
              case 'oil water'
                fluidModel = OilWaterFixedStressFluidModel(model.G, rock, fluid, varargin{:});
              case 'blackoil'
                fluidModel = BlackOilFixedStressFluidModel(model.G, rock, fluid, varargin{:});
              otherwise
                error('fluidModelType not recognized.');
            end
            % Validate model for simulation
            fluidModel = fluidModel.validateModel();
        end



        function [state, report] = stepFunction(model, state, state0, dt, ...
                                                drivingForces, linsolve, ...
                                                nonlinsolve, iteration, ...
                                                varargin)

            fluidModel = model.fluidModel;
            mechModel = model.mechModel;

            % The state variable has two parts, the mechanical part mstate
            % and the fluid part wstate, see the synchronization member
            % functions in the MechFluidSplitModel.

            mstate0 = model.syncMStateFromState(state0);
            wstate0 = model.syncWStateFromState(state0);

            % We get the fluid pressure which is an input for the mechanical system
            fluidp = fluidModel.getProp(state, 'pressure');

            mechsolver = model.mech_solver;

            if model.verbose
                fprintf('=== Splitting scheme: step %d\n', iteration);
            end

            % We solve the mechanical system
            [mstate, mreport] = mechsolver.solveTimestep(mstate0, dt, mechModel, ...
                                                          'fluidp', fluidp);

            if model.verbose
                fprintf('Solved for mechanics (%d iterations)\n', mreport.Iterations);
            end

            wdrivingForces = drivingForces; % The main model gets the well controls

            % We update the global (mech+fluid) state with the newly computed mechanical
            % state
            state = model.syncStateFromMState(state, mstate);

            % We compute the change of effective volume induced by the mechanical
            % part. These are integrated in the driving force for the fluid
            % system.
            wdrivingForces.fixedStressTerms.new = computeMechTerm(model, state);
            wdrivingForces.fixedStressTerms.old = computeMechTerm(model, state0);

            forceArg = fluidModel.getDrivingForces(wdrivingForces);

            fluidsolver = model.fluid_solver;

            % We solve the fluid system.
            [wstate, wreport] = fluidsolver.solveTimestep(wstate0, dt, fluidModel, ...
                                                          forceArg{:});

            if model.verbose
                fprintf('Solved for fluid (%d iterations)\n', wreport.Iterations);
            end

            % We update the global state with the newly computed fluid state
            state = model.syncStateFromWState(state, wstate);

            % We check the convergence by looking back that the mechanical
            % system and check the residual for these equations.
            fluidp = fluidModel.getProp(state, 'pressure');
            drivingForces.fluidp = fluidp;
            mstate = model.syncMStateFromState(state);
            problem = mechModel.getEquations(mstate0, mstate, dt, drivingForces, ...
                                                      'resOnly', true, ...
                                                      'iteration', iteration);
            [convergence, values, names] = model.checkConvergence(problem);
            failureMsg = '';
            failure = false;
            isConverged = all(convergence);

            if model.verbose
                fprintf(['Value of residual for mechanical part: ' ...
                         '%g\n'], values);
                if convergence
                    fprintf('Convergence reached\n');
                else
                    fprintf(['Convergence not reached (tolerance value = ' ...
                             '%g)\n'], model.splittingTolerance);
                end
            end

            report = mechModel.makeStepReport('Failure',      failure, ...
                                              'FailureMsg',   failureMsg, ...
                                              'Converged',    isConverged, ...
                                              'Residuals',    values, ...
                                              'ResidualsConverged', ...
                                              convergence);
        end

        function fixedStressTerms = computeMechTerm(model, state)
        % Compute the change in pore volume induced by the mechanical system. In a fixed
        % stress splitting method, the relation between volume change and
        % pressure is obtained by assuming that the total stress is
        % preserved.

            mechModel = model.mechModel;
            %stress = mechModel.operators.stress*(state.xd);
            
            uvec = state.u;
            if isa( state.xd, 'ADI')
               uvec = double2ADI(uvec, state.xd);
            end
            uvec(~mechModel.operators.isdirdofs) = state.xd;
            stress = mechModel.operators.global_stress * uvec;
            
            p = state.pressure;

            invCi = mechModel.mech.invCi;
            % invCi is the tensor equal to $C^{-1}I$ where $I$ is the identity tensor.
            griddim = model.G.griddim;

            pTerm = mechModel.operators.trace(invCi); % could have been
                                                            % computed and
                                                            % stored...

            if griddim == 3
                pI = bsxfun(@times, p, [1, 1, 1, 0, 0, 0]);
                nlin = 6;
            else
                pI = bsxfun(@times, p, [1, 1, 0]);
                nlin = 3;
            end
            stress = reshape(stress, nlin, []);
            stress = stress';
            totalStress = stress - pI;

            sTerm = sum(invCi.*totalStress, 2);

            fixedStressTerms.pTerm = pTerm; % Compressibility due to mechanics
            fixedStressTerms.sTerm = sTerm; % Volume change due to mechanics

            % Short explanation about the terms computed above:
            %
            % We have $\sigma_T = C\epsi - p I$ where $\sigma_T$ is total
            % stress
            % Hence, $tr(\epsi) = tr(C^{-1}\sigma_T) + p tr(C^{-1})$ and the
            % value sTerm and pTerm corresponds to $tr(\epsi) =
            % tr(C^{-1}\sigma_T)$ and $tr(C^{-1})$ respectively.
            %
            % Note that $tr(C^{-1}\sigma_T) = I:(C^{-1}\sigma_T) = invCi:\sigma_T$
        end


    end
end

classdef MechFluidFixedStressSplitModel < MechFluidSplitModel


    methods
        function model = MechFluidFixedStressSplitModel(G, rock, fluid, mech_problem, varargin)
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
        end


        function [state, report] = stepFunction(model, state, state0, dt, ...
                                                drivingForces, linsolve, ...
                                                nonlinsolve, iteration, ...
                                                varargin)

            fluidModel = model.fluidModel;
            mechModel = model.mechModel;

            % Solve the mechanic equations
            mstate0 = model.syncMStateFromState(state0);
            wstate0 = model.syncWStateFromState(state0);
            state = state0;

            fluidp = fluidModel.getProp(wstate0, 'pressure');
            mechsolver = model.mech_solver;

            [mstate, mreport] = mechsolver.solveTimestep(mstate0, dt, mechModel, ...
                                                          'fluidp', fluidp);

            % Solve the fluid equations
            wdrivingForces = drivingForces; % The main model gets the well controls

            state = model.syncStateFromMState(state, mstate);

            wdrivingForces.fixedStressTerms.new = computeMechTerm(model, state);
            wdrivingForces.fixedStressTerms.old = computeMechTerm(model, state0);

            forceArg = fluidModel.getDrivingForces(wdrivingForces);

            fluidsolver = model.fluid_solver;
            [wstate, wreport] = fluidsolver.solveTimestep(wstate0, dt, fluidModel, ...
                                                          forceArg{:});

            state = model.syncStateFromMState(state, mstate);
            state = model.syncStateFromWState(state, wstate);

            report = model.makeStepReport( 'Failure',         false, 'Converged', ...
                                           true, 'FailureMsg',      '', ...
                                           'Residuals',       0 );

        end


        function fixedStressTerms = computeMechTerm(model, state)
            stress = state.stress;
            p = state.pressure;

            invCi = model.mechModel.mech.invCi;
            griddim = model.G.griddim;

            pTerm = sum(invCi(:, 1 : griddim), 2); % could have been computed and stored...

            if griddim == 3
                cvoigt = [1, 1, 1, 0.5, 0.5, 0.5];
            else
                cvoigt = [1, 1, 0.5];
            end
            stress = bsxfun(@times, stress, cvoigt);
            sTerm = sum(invCi.*stress, 2);

            fixedStressTerms.pTerm = pTerm; % Compressibility due to mechanics
            fixedStressTerms.sTerm = sTerm; % Volume change due to mechanics

        end

    end
end

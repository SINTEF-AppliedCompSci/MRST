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

            if iteration == 1
                state_prev = state0;
            else
                state_prev = state;
            end

            fluidp = fluidModel.getProp(wstate0, 'pressure');
            mechsolver = model.mech_solver;

            [mstate, mreport] = mechsolver.solveTimestep(mstate0, dt, mechModel, ...
                                                          'fluidp', fluidp);

            % Solve the fluid equations
            wdrivingForces = drivingForces; % The main model gets the well controls

            state = model.syncStateFromMState(state, mstate);

            wdrivingForces.fixedStressTerms.new = computeMechTerm(model, state);
            wdrivingForces.fixedStressTerms.old = computeMechTerm(model, state_prev);

            forceArg = fluidModel.getDrivingForces(wdrivingForces);

            fluidsolver = model.fluid_solver;
            [wstate, wreport] = fluidsolver.solveTimestep(wstate0, dt, fluidModel, ...
                                                          forceArg{:});

            state = model.syncStateFromMState(state, mstate);
            state = model.syncStateFromWState(state, wstate);

            report.Converged = false;
            report.Failure   = false;
            report.Residuals = [];

            if iteration > 1
                [incAbs, incVarNames] = model.computeNormIncrements(state_prev, ...
                                                                  state);
                if norm(incAbs, inf) < model.nonlinearTolerance
                    report.Converged = true;
                end
            end

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

        function [incAbs, incVarNames] = computeNormIncrements(model, state_prev, state)

            mechfds  = model.mechfds;
            fluidfds = model.fluidfds;
            fluidfds = model.stripVars(fluidfds, 'wellSol');

            facilitymodel =  model.FacilityModel;
            wellSol = state.wellSol;
            wellSol_prev = state_prev.wellSol;
            wellVars = facilitymodel.getPrimaryVariableNames();
            nwellVars = numel(wellVars);
            actIx = facilitymodel.getIndicesOfActiveWells();
            nW = numel(actIx);

            nInc = numel(mechfds) + numel(fluidfds) + nwellVars;
            incAbs = zeros(nInc, 1);
            incVarNames = cell(1, nInc);

            varNo = 1;

            for varMechNo = 1 : numel(mechfds)
                incVarNames{varNo} = mechfds{varMechNo};
                incAbs(varNo) = norm(model.mechModel.getProp(state, ...
                                                             mechfds{varMechNo}) ...
                                     - model.mechModel.getProp(state_prev, ...
                                                               mechfds{varMechNo}), ...
                                     Inf);
                varNo = varNo + 1;
            end

            for varFluidNo = 1 : numel(fluidfds)
                incVarNames{varNo} = fluidfds{varFluidNo};
                incAbs(varNo) = norm(model.fluidModel.getProp(state, ...
                                                              fluidfds{varFluidNo}) ...
                                     - model.fluidModel.getProp(state_prev, ...
                                                                fluidfds{varFluidNo}), ...
                                     Inf);
                varNo = varNo + 1;
            end

            for varWellNo = 1 : nwellVars
                wf = wellVars{varWellNo};
                incVarNames{varNo} = wf;
                isVarWell = facilitymodel.getWellVariableMap(wf);
                for wNo = 1 : nW
                    subs = (isVarWell == wNo);
                    if any(subs)
                        incAbs(varNo) = incAbs(varNo) + abs ...
                            (facilitymodel.WellModels{wNo}.getProps(wellSol(wNo), ...
                                                                    wf)- ...
                             facilitymodel.WellModels{wNo}.getProps(wellSol_prev(wNo), ...
                                                                    wf));
                    end
                end
                varNo = varNo + 1;
            end


        end


    end
end

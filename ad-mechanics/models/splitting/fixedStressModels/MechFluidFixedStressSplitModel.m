classdef MechFluidFixedStressSplitModel < MechFluidSplitModel
% Base class to implement fixed-stress splitting.
%
% Different fluid model can be used. Each of these will use a derived class
% from this one.

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


            % problemFluid = fluidModel.getEquations(state0, state, dt, ...
            %                                                wdrivingForces, ...
            %                                                'iteration', ...
            %                                                iteration, 'resOnly', ...
            %                                                true);
            
            fcModel = model.fullyCoupledModel;
            problem = fcModel.getEquations(state0, state, dt, drivingForces, ...
                                                   'resOnly', true, 'iteration', ...
                                                   iteration);
            
            [convergence, values, names] = fcModel.checkConvergence(problem);            
            failureMsg = '';
            failure = false;
            isConverged = all(convergence);

            report = fcModel.makeStepReport( 'Failure',      failure, ...
                                             'FailureMsg',   failureMsg, ...
                                             'Converged',    isConverged, ...
                                             'Residuals',    values, ...
                                             'ResidualsConverged', ...
                                             convergence);

        end

        function [problem] = getFullyCoupledEquations()
        % The residual of the fully coupled equations are computed to test convergence.
            error('Base class function not meant for direct use.');
        end

        function fixedStressTerms = computeMechTerm(model, state)
            stress = state.stress;
            p = state.pressure;

            invCi = model.mechModel.mech.invCi;
            % invCi is the tensor equal to $C^{-1}I$ where $I$ is the identity tensor.
            griddim = model.G.griddim;

            pTerm = model.mechModel.operators.trace(invCi); % could have been
                                                            % computed and
                                                            % stored...

            if griddim == 3
                cvoigt = [1, 1, 1, 0.5, 0.5, 0.5];
                pI = bsxfun(@times, p, [1, 1, 1, 0, 0, 0]);
            else
                cvoigt = [1, 1, 0.5];
                pI = bsxfun(@times, p, [1, 1, 0]);
            end
            totalStress = stress - pI;
            totalStress = bsxfun(@times, totalStress, cvoigt);
            
            sTerm = sum(invCi.*totalStress, 2); 

            fixedStressTerms.pTerm = pTerm; % Compressibility due to mechanics
            fixedStressTerms.sTerm = sTerm; % Volume change due to mechanics

            % We have $\sigma_T = C\epsi - p I$ where $\sigma_T$ is total
            % stress
            % Hence, $tr(\epsi) = tr(C^{-1}\sigma_T) + p tr(C^{-1})$ and the
            % value sTerm and pTerm corresponds to $tr(\epsi) =
            % tr(C^{-1}\sigma_T)$ and $tr(C^{-1})$ respectively.
            %
            % Note that $tr(C^{-1}\sigma_T) = I:(C^{-1}\sigma_T) = invCi:\sigma_T$
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

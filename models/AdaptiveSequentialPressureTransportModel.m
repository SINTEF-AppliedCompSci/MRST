classdef AdaptiveSequentialPressureTransportModel < SequentialPressureTransportModel
    % Sequential meta-model which solves pressure and transport using a fixed
    % flux splitting
    properties
        GF
        G0
        rockF
    end
    
    methods
        function model = AdaptiveSequentialPressureTransportModel(pressureModel, transportModel, GF, G0, rockF, varargin)
            model = model@SequentialPressureTransportModel(pressureModel, transportModel, varargin{:});
            model.GF = GF;
            model.G0 = G0;
            model.rockF = rockF;
        end
        
        
        function [state, report] = stepFunction(model, state, state0, dt,...
                                                drivingForces, linsolve, nonlinsolve,...
                                                iteration, varargin)
            
            if isfield(state, 'model')
                model = state.model;
                drivingForces = state.forces;
            end
            [state, pressureReport, transportReport, pressure_ok, transport_ok, forceArg] =...
                model.solvePressureTransport(state, state0, dt, drivingForces, iteration);

            converged = pressure_ok && transport_ok;
            if converged && ~model.stepFunctionIsLinear
                [converged, values, state] = checkOuterConvergence(model, state, state0, dt, drivingForces, iteration);
            else
                % Need to have some value in the report
                values = pressureReport.StepReports{end}.NonlinearReport{end}.Residuals(1);
            end
            failure = false;
            FailureMsg = '';
            if ~pressure_ok
                converged = converged && false;
            end
            report = model.makeStepReport(...
                                    'Failure',         failure, ...
                                    'Converged',       all(converged), ...
                                    'FailureMsg',      FailureMsg, ...
                                    'ResidualsConverged', converged, ...
                                    'Residuals',       values ...
                                    );
                                
            report.PressureSolver =  pressureReport;
            report.TransportSolver = transportReport;
            
            if model.reupdatePressure && converged
                state = ...
                    psolver.solveTimestep(state0, dt, model.pressureModel,...
                            'initialGuess', state, ...
                            forceArg{:});
                [~, state] = model.transportModel.getEquations(state0, state, dt, drivingForces, 'resOnly', true, 'iteration', inf);
            end
        end
        
        function [state, pressureReport, transportReport, pressure_ok, transport_ok, forceArg, model] = solvePressureTransport(model, state, state0, dt, drivingForces, iteration)
           % Solve pressure and transport sequentially
            psolver = model.pressureNonLinearSolver;
            tsolver = model.transportNonLinearSolver;
            if iteration > 1
                for i = 1:numel(model.pressureModel.FacilityModel.WellModels)
                    model.pressureModel.FacilityModel.WellModels{i}.doUpdatePressureDrop = false;
                end
            end
            
            residual = 0;
            if isfield(state.wellSol, 'flux')
                problem = model.transportModel.getEquations(state0, state, ...
                                           dt, drivingForces, 'resOnly', true);
                for eqNo = 1:numel(problem)
                    residual = residual + abs(problem.equations{eqNo});
                end 
            end

            G = model.transportModel.G;
            wc = [drivingForces.W.cells];
            residual(wc) = 0;
            if isfield(G, 'refined') && 0
                residual(G.refined) = false;
            end

            tol = 1e-2;
            cells = find(residual > tol);
            if any(cells)
                [model, state, state0, drivingForces] = model.refineModel(cells, state, state0, drivingForces);
            end
            state.G = model.transportModel.G;
            state.model = model;
            state.forces = drivingForces;
            
            % Get the forces used in the step
            forceArg = model.pressureModel.getDrivingForces(drivingForces);
            
            % First, solve the pressure using the pressure nonlinear
            % solver.
            [state, pressureReport] = ...
                psolver.solveTimestep(state0, dt, model.pressureModel,...
                            'initialGuess', state, ...
                            forceArg{:});
            pressure_ok = pressureReport.Converged || psolver.continueOnFailure;
            
            if pressure_ok

                if ~isempty(drivingForces.bc)
                    isDir = strcmpi(drivingForces.bc.type, 'pressure');
                    if any(isDir)
                        % Convert Dirichlet boundary conditions to flux
                        % boundary conditions for the transport
                        transportForces = drivingForces;
                        G = model.pressureModel.G;
                        dirFace = transportForces.bc.face(isDir);
                        q = sum(state.flux(dirFace, :), 2);
                        sgn = 1 - 2*(G.faces.neighbors(dirFace, 2) == 0);
                        transportForces.bc.value(isDir) = sgn.*q;
                        [transportForces.bc.type{isDir}] = deal('resflux');
                        forceArg = model.transportModel.getDrivingForces(transportForces);
                    end
                end
                
                state.timestep = dt;
                state.pressure_full = state.pressure;
                
                % If pressure converged, we proceed to solve the transport
                [state, transportReport] = ...
                    tsolver.solveTimestep(state0, dt, model.transportModel,...
                                'initialGuess', state, ...
                                forceArg{:});
                transport_ok = transportReport.Converged;
                
                % Map quantities from 
                
            else
                transport_ok = false;
                transportReport = [];
            end 
        end
        
        function [model, state, state0, forces] = refineModel(model, cells, state, state0, forces)
            
            tm = model.transportModel;
            G  = tm.G;
            
             % Refine grid
            [G, mappings] = refineGrid(G, model.G0, model.GF, cells);
            
            sF  = state.s(mappings.fine2old,:);
            bF  = state.bfactor(mappings.fine2old,:);
            pvF = poreVolume(model.GF, model.rockF);
            
            sF0 = state0.s(mappings.fine2old,:);
            bF0 = state0.bfactor(mappings.fine2old,:);
            
            S = sparse(mappings.fine2new, (1:model.GF.cells.num)', 1);
            state.s  = S*(sF.*bF.*pvF)./(S*(bF.*pvF));
            state0.s = S*(sF0.*bF0.*pvF)./(S*(bF0.*pvF));
            
            pF             = state.pressure(mappings.fine2old);
            state.pressure = S*pF./sum(S>0,2);
            
            pF0             = state0.pressure(mappings.fine2old);
            state0.pressure = S*pF0./sum(S>0,2);
            
            rock = makeRock(G, 1, 1);
            rock.perm = model.rockF.perm(mappings.new2fine,:);
            rock.poro = model.rockF.poro(mappings.new2fine,:);
            
            op = setupOperatorsTPFA(G, rock);
            
            model.transportModel.G = G;
            model.transportModel.rock = rock;
            model.transportModel.operators = op;
            model.pressureModel.G = G;
            model.pressureModel.rock = rock;
            model.pressureModel.operators = op;
            
            state.cells = cells;

            if ~isempty(forces.W)
                W = forces.W;
                for wNo = 1:numel(W)
                    W(wNo).cells = mappings.old2new(W(wNo).cells);
%                     W(wNo).cells = mappings.new2
                end
                forces.W = W;
            end
            
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

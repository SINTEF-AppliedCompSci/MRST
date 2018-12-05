classdef AdaptiveSequentialPressureTransportModel < SequentialPressureTransportModel
    % Sequential meta-model which solves pressure and transport using a fixed
    % flux splitting
    properties
        GF
        G0
        rockF
        coarseTransportModel
        fineTransportModel
        plotProgress
    end
    
    methods
        function model = AdaptiveSequentialPressureTransportModel(pressureModel, transportModel, GF, G0, rockF, varargin)
            model = model@SequentialPressureTransportModel(pressureModel, transportModel, varargin{:});
            model.GF = GF;
            model.G0 = G0;
            model.rockF = rockF;
            model.fineTransportModel = transportModel;
            model.coarseTransportModel = upscaleModelTPFA(transportModel, G0.partition);
            model.plotProgress = true;
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
             % Get the forces used in the step
            forceArg = model.pressureModel.getDrivingForces(drivingForces);
            
            % First, solve the pressure using the pressure nonlinear
            % solver.
            [state, pressureReport] = ...
                psolver.solveTimestep(state0, dt, model.pressureModel,...
                            'initialGuess', state, ...
                            forceArg{:});
            pressure_ok = pressureReport.Converged || psolver.continueOnFailure;
            statePressure = state;
            
            if pressure_ok
                % Assemble fine-scale transport residual
                
                stateUpsc  = upscaleState(model.coarseTransportModel, model.fineTransportModel, state);
                stateUpsc0 = upscaleState(model.coarseTransportModel, model.fineTransportModel, state0);
                drivingForcesUpsc = model.mapForces(drivingForces, model.G0);
                
                residual = 0;
                if isfield(state.wellSol, 'flux')
                    problem = model.coarseTransportModel.getEquations(stateUpsc0, stateUpsc, ...
                                               dt, drivingForcesUpsc, 'resOnly', true);
                    for eqNo = 1:numel(problem)
                        residual = residual + abs(problem.equations{eqNo});
                    end 
                end

                % Make sure we dont refine well cells (they should already
                % be at a suitable refinemnet level)
                wc = model.G0.partition([drivingForces.W.cells]);
                G = model.G;
                residual(wc) = 0;
                if isfield(G, 'refined') && 0
                    residual(G.refined) = 0;
                end
            
            
                tol = 1e-2;
                cells = abs(residual) > tol;
                if 1
                [model, state, state0, drivingForces, mappings] = model.refineModel(cells, state, state0, drivingForces);
                end
                transportForces = drivingForces;
                state.G = model.transportModel.G;
                state.forces = transportForces;

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
                    end
                end
                forceArg = model.transportModel.getDrivingForces(transportForces);
                
                state.timestep = dt;
                state.pressure_full = state.pressure;
                
                % If pressure converged, we proceed to solve the transport
                [state, transportReport] = ...
                    tsolver.solveTimestep(state0, dt, model.transportModel,...
                                'initialGuess', state, ...
                                forceArg{:});
                transport_ok = transportReport.Converged;
                
                if 1
                stateFine = statePressure;
                stateFine.s = state.s(mappings.fine2new,:);
                stateFine.G = state.G;
                state = stateFine;
                end
                
            else
                transport_ok = false;
                transportReport = [];
            end 
        end
        
        function [model, state, state0, forces, mappings] = refineModel(model, cells, state, state0, forces)
            
            tm = model.transportModel;
            G  = tm.G;
            
             % Refine grid
            [G, mappings, partition] = refineGrid(model.G0, model.G0, model.GF, cells);
            
            model.transportModel = upscaleModelTPFA(model.transportModel, partition);
            
%             rock = makeRock(G, 1, 1);
%             rock.perm = model.rockF.perm(mappings.new2fine,:);
%             rock.poro = model.rockF.poro(mappings.new2fine,:);
%             op = setupOperatorsTPFA(G, rock);
%             model.transportModel.G = G;
%             model.transportModel.rock = rock;
%             model.transportModel.operators = op;
            

%             pvC = model.transportModel.operators.pv;
%             pvF = model.fineTransportModel.operators.pv;
%             
%             S = sparse(mappings.fine2new, (1:model.GF.cells.num)', 1);
%             
%             state.pressure  = S*(pvF.*state.pressure)./pvC;
%             state0.pressure = S*(pvF.*state0.pressure)./pvC;
%             bF0 = state0.bfactor;
%             bF = state.bfactor;
%             bC = S*(pvF.*state.bfactor)./pvC;
%             bC0 = S*(pvF.*state0.bfactor)./pvC;
% 
%             if isfield(state, 'flux')
%                 cfsign = fineToCoarseSign(G);
%                 cfacesno = rldecode(1:G.faces.num, diff(G.faces.connPos), 2) .';
%                 newflux = zeros(G.faces.num, size(state.flux, 2));
%                 for i = 1:size(state.flux, 2)
%                     newflux(:, i)   = accumarray(cfacesno, state.flux(G.faces.fconn, i) .* cfsign);
%                 end
%                 state.flux = newflux;
%             end
%             
% %             state.s  = S*(sF.*bF.*pvF)./(S*(bF.*pvF));
% %             state0.s = S*(sF0.*bF0.*pvF)./(S*(bF0.*pvF));
%             state.s = S*bsxfun(@times, pvF, bF.*state.s)./bsxfun(@times, pvC, bC);
%             
%             state.bfactor = bC;
%             
%             state0.s = S*bsxfun(@times, pvF, bF0.*state0.s)./bsxfun(@times, pvC, bC0);
%             state0.bfactor = bC0;
            
%             [pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.transportModel.fluid, p, p0);

            state = upscaleState(model.transportModel, model.fineTransportModel, state);
            state.G = G;
            state0 = upscaleState(model.transportModel, model.fineTransportModel, state0);
%             sF  = state.s(mappings.fine2old,:);
%             bF  = state.bfactor(mappings.fine2old,:);
%             pvF = poreVolume(model.GF, model.rockF);
%             
%             sF0 = state0.s(mappings.fine2old,:);
%             bF0 = state0.bfactor(mappings.fine2old,:);
%             
%             S = sparse(mappings.fine2new, (1:model.GF.cells.num)', 1);
%             state.s  = S*(sF.*bF.*pvF)./(S*(bF.*pvF));
%             state0.s = S*(sF0.*bF0.*pvF)./(S*(bF0.*pvF));
%             
%             pF             = state.pressure(mappings.fine2old);
%             state.pressure = S*pF./sum(S>0,2);
%             
%             pF0             = state0.pressure(mappings.fine2old);
%             state0.pressure = S*pF0./sum(S>0,2);
%             
%             rock = makeRock(G, 1, 1);
%             rock.perm = model.rockF.perm(mappings.new2fine,:);
%             rock.poro = model.rockF.poro(mappings.new2fine,:);
%             
%             op = setupOperatorsTPFA(G, rock);
%             
%             model.transportModel.G = G;
%             model.transportModel.rock = rock;
%             model.transportModel.operators = op;
%             model.pressureModel.G = G;
%             model.pressureModel.rock = rock;
%             model.pressureModel.operators = op;
%             
%             state.cells = cells;


            forces = model.mapForces(forces, G);
%             if ~isempty(forces.W)
%                 W = forces.W;
%                 for wNo = 1:numel(W)
%                     W(wNo).cells = mappings.fine2new(W(wNo).cells);
% %                     W(wNo).cells = mappings.new2
%                 end
%                 forces.W = W;
%             end
            
            if model.plotProgress
                figure(1)
                clf;
                plotCellData(G, state.s(:,1), 'edgecolor', 'none');
                plotGrid(G, 'facec', 'none', 'edgealpha', 0.3);
                colormap(jet)
                axis equal tight
                caxis([0.2,0.8])
                drawnow()
            end
            
        end
        
        function forces = mapForces(model, forces, G)
            if ~isempty(forces.W)
                W = forces.W;
                for wNo = 1:numel(W)
                    W(wNo).cells = G.partition(W(wNo).cells);
                end
                forces.W = W;
            end
            
        end
        
        function state = upscaleState(model);
            
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

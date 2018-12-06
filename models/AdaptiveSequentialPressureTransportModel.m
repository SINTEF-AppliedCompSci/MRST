classdef AdaptiveSequentialPressureTransportModel < SequentialPressureTransportModel
    % Sequential meta-model which solves pressure and transport using a fixed
    % flux splitting
    properties
        G0
        coarseTransportModel
        fineTransportModel
        plotProgress
        storeGrids
    end
    
    methods
        function model = AdaptiveSequentialPressureTransportModel(pressureModel, transportModel, G0, varargin)
            model = model@SequentialPressureTransportModel(pressureModel, transportModel, varargin{:});
            model.G0 = G0;
            model.fineTransportModel = transportModel;
            model.coarseTransportModel = upscaleModelTPFA(transportModel, G0.partition);
            model.plotProgress = false;
            model.storeGrids   = false;
            model = merge_options(model, varargin{:});
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
%                     for eqNo = 1:numel(problem)
%                         residual = residual + abs(problem.equations{eqNo});
%                     end 
                    residual = horzcat(problem.equations{:});
                end

                % Make sure we dont refine well cells (they should already
                % be at a suitable refinemnet level)
                wc = model.coarseTransportModel.G.partition([drivingForces.W.cells]);
                residual(wc) = 0;
            
                tol = 1e-2;
                cells = abs(residual) > tol;
                [model, state, state0, drivingForces] ...
                    = model.refineTransportModel(cells, state, state0, drivingForces);
                transportForces = drivingForces;
                state.forces    = transportForces;

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
                
                % Map transport quantities onto fine grid to prepare for
                % next pressure solve
                varNames = getTransportVarNames();
                p     = model.transportModel.G.partition;
                st    = state;
                state = statePressure;
                state.G = st.G;
                for vNo = 1:numel(varNames)
                    vn = varNames{vNo};
                    if isfield(state, vn)
                        state.(vn) = st.(vn)(p,:);
                    end
                end
                
            else
                transport_ok = false;
                transportReport = [];
            end 
        end
        
        function [model, state, state0, forces, mappings] = refineTransportModel(model, cells, state, state0, forces)
            
            GF = model.fineTransportModel.G;
            GC = model.coarseTransportModel.G;
            
             % Refine grid
            [G, mappings, partition] = refineGrid(GC, GC, GF, cells);
            model.transportModel = upscaleModelTPFA(model.transportModel, partition);
            model.transportModel.G.cells.ghost = false(G.cells.num,1);
            
            state  = model.upscaleState(state);
            state0 = model.upscaleState(state0);
            
            if model.storeGrids
                state.G = G;
            end

            forces = model.mapForces(forces, G);
            
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
        
        function forces = mapForces(model, forces, G) %#ok
            if ~isempty(forces.W)
                W = forces.W;
                for wNo = 1:numel(W)
                    W(wNo).cells = G.partition(W(wNo).cells);
                end
                forces.W = W;
            end
            
        end
        
        function state = upscaleState(model, state)
            
            % Fine and coarsened grids
            GF = model.fineTransportModel.G;
            G  = model.transportModel.G;
            
            % Pore volumes
            pvf = model.fineTransportModel.operators.pv;
            pv  = model.transportModel.operators.pv;
            
            % Summation matrix
            S = sparse(G.partition, (1:GF.cells.num)', 1);
            
            % 
            pvbf           = bsxfun(@times, pvf, state.bfactor);
            state.bfactor  = S*pvbf./pv;
            state.pressure = S*(pvf.*state.pressure)./pv;
            state.s        = S*(pvbf.*state.s)./(S*pvbf);

            cFsign   = fineToCoarseSign(G);
            cFacesno = rldecode(1:G.faces.num, diff(G.faces.connPos), 2) .';
            newFlux = zeros(G.faces.num, size(state.flux, 2));
            for i = 1:size(state.flux, 2)
                newFlux(:, i)   = accumarray(cFacesno, state.flux(G.faces.fconn, i) .* cFsign);
            end
            state.flux = newFlux;
        end

    end
end

function varNames = getTransportVarNames()

    varNames = {'s'};
    
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

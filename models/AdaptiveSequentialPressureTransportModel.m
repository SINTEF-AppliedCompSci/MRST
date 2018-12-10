classdef AdaptiveSequentialPressureTransportModel < SequentialPressureTransportModel
    % Sequential meta-model which solves pressure and transport using a fixed
    % flux splitting
    properties
        G0
        coarseTransportModel
        fineTransportModel
        plotProgress
        storeGrids
        isDG
    end
    
    methods
        function model = AdaptiveSequentialPressureTransportModel(pressureModel, transportModel, G0, varargin)
            model = model@SequentialPressureTransportModel(pressureModel, transportModel, varargin{:});
            model.G0 = G0;
            model.fineTransportModel = transportModel;
            model.coarseTransportModel = model.upscaleTransportModelTPFA(G0.partition);
            model.transportModel = model.coarseTransportModel;
            model.coarseTransportModel.G.cells.refined = G0.cells.refined;
            model.coarseTransportModel.G.oldPartition = G0.partition;
            
            model.transportModel = model.coarseTransportModel;
            model.plotProgress = false;
            model.storeGrids   = false;
            model.isDG = isa(model.transportModel, 'TransportBlackOilModelDG');
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
            [pressureState, pressureReport] = ...
                psolver.solveTimestep(state0, dt, model.pressureModel,...
                            'initialGuess', state, ...
                            forceArg{:});
            pressure_ok = pressureReport.Converged || psolver.continueOnFailure;
            
            if pressure_ok
                % Assemble transport residual
                transportState    = model.upscaleState(pressureState);
                transportState0   = model.upscaleState(state0);
                drivingForcesUpsc = model.mapForces(drivingForces, model.G0);
                
                problem = model.transportModel.getEquations(transportState0, transportState, ...
                                               dt, drivingForcesUpsc, 'resOnly', true);
                residual = horzcat(problem.equations{:});
                if isprop(model.transportModel, 'disc')
                    ix = rldecode((1:model.transportModel.G.cells.num)', transportState.nDof, 1);
                    residual = accumarray(ix, abs(residual));
                end

                % Make sure we dont refine already refined cells
                if isfield(model.transportModel.G.cells, 'refined')
                    residual(model.transportModel.G.cells.refined) = 0;
                end
            
                tol = 1e-2;
                cells = find(abs(residual) > tol);
                [model, transportState, transportState0, drivingForces] ...
                    = model.refineTransportModel(cells, pressureState, state0, drivingForces);
                transportForces = drivingForces;

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
                
                transportState.timestep = dt;
                transportState.pressure_full = transportState.pressure;
                
                % If pressure converged, we proceed to solve the transport
                [transportState, transportReport] = ...
                    tsolver.solveTimestep(transportState0, dt, model.transportModel,...
                                'initialGuess', transportState, ...
                                forceArg{:});
                transport_ok = transportReport.Converged;
                
                % Map transport quantities onto fine grid to prepare for
                % next pressure solve
                state = model.refineState(transportState, pressureState);
                
            else
                transport_ok = false;
                transportReport = [];
            end 
        end
        
        function [model, transportState, transportState0, forces] ...
                = refineTransportModel(model, cells, transportState, transportState0, forces)
            
            GF = model.fineTransportModel.G;
            GC = model.coarseTransportModel.G;
            if isfield(transportState0, 'G') && isfield(transportState0.G, 'parent')
                G = transportState0.G;
            else
                G = GC;
            end
            
            % Refine grid
            mappings = getRefinementMappings(G, GC, GF, cells);
            tm       = model.upscaleTransportModelTPFA(mappings.newPartition);
            
            model.transportModel = tm;
            transportState.pv = model.transportModel.operators.pv;
            transportState  = model.upscaleState(transportState);
            transportState0 = model.upscaleState(transportState0);
            
            if model.storeGrids
                transportState.G = tm.G;
            end

            forces = model.mapForces(forces, tm.G);
            
            if model.plotProgress
                figure(1)
                clf;
                plotCellData(tm.G, transportState.s(:,1), 'edgecolor', 'none');
                plotGrid(tm.G, 'facec', 'none', 'edgealpha', 0.3);
                colormap(jet)
                axis equal tight
%                 caxis([0.2,0.8])
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
            
            % Summation matrix
            S = S*sparse((1:GF.cells.num)', state.transportState.G.partition, 1);
%             SS = S*SS';
            
            pvPrev = state.transportState.pv;
            pvb = bsxfun(@times, pv, state.bfactor);
            state.s = S*(pvPrev.*state.transportState.s)./(S*pvPrev);
            
            ts = state.transportState;
            
            if isa(model.transportModel, 'TransportBlackOilModelDG')
                disc = model.transportModel.disc;
                state.degree = repmat(disc.degree, G.cells.num, 1);
                state = disc.updateDofPos(state);
                sdof = zeros(sum(state.nDof), size(state.s,2));
                for dofNo = 1:numel(disc.basis.nDof)
                    ix     = disc.getDofIx(ts, dofNo, 1:size(ts.s,1));
                    ixNew  = disc.getDofIx(state, dofNo, 1:G.cells.num);
                    keepRows = state.nDof >= dofNo;
                    keepCols = ts.nDof >= dofNo;
                    Stmp = S(keepRows, keepCols);
                    
                    sdof(ixNew,:) = Stmp*(pvPrev.*ts.sdof(ix,:))./(Stmp*pvPrev);
                end
                state.sdof = sdof;
                state.degree = repmat(model.transportModel.disc.degree, G.cells.num, 1);            
            	state        = model.transportModel.disc.assignDofFromState(state);
            end
            
%             if 
%             state.degree   = repmat(model.transportModel.disc.degree, G.cells.num, 1);

            cFsign   = fineToCoarseSign(G);
            cFacesno = rldecode(1:G.faces.num, diff(G.faces.connPos), 2).';
            newFlux = zeros(G.faces.num, size(state.flux, 2));
            for i = 1:size(state.flux, 2)
                newFlux(:, i)   = accumarray(cFacesno, state.flux(G.faces.fconn, i) .* cFsign);
            end
            state.flux = newFlux;
        end
        
        function transportModel = upscaleTransportModelTPFA(model, partition, varargin)
            
            transportModel = upscaleModelTPFA(model.fineTransportModel, partition, varargin{:});
            G = coarsenCellDimensions(transportModel.G);
            G.cells.ghost   = false(G.cells.num, 1);
            G.cells.refined = accumarray(G.partition, ones(G.parent.cells.num,1)) == 1;
            transportModel.G = G;
            
            if model.isDG
                [jt, ot, mt] = deal(Inf);
                transportModel.disc = DGDiscretization(transportModel     , ...
                                    'degree'               , model.transportModel.disc.degree, ...
                                    'basis'                , 'legendre'   , ...
                                    'useUnstructCubature'  , false        , ...
                                    'jumpTolerance'        , jt           , ...
                                    'outTolerance'         , ot           , ...
                                    'outLimiter'           , 'kill'       , ...
                                    'meanTolerance'        , mt           , ...
                                    'limitAfterConvergence', false        , ...
                                    'plotLimiterProgress'  , false        );
            end
            
        end
        
        function state = refineState(model, transportState, pressureState)
            
            [varNames, dofVarNames] = getTransportVarNames();
            p     = model.transportModel.G.partition;
            st    = transportState;
            
            
            state = pressureState;
%             state.wellSol = transportState.wellSol;
            
%             st.pvPrev = model.transportModel.operators.pv;
%             st.oldPartition = model.transportModel.G.partition;

            state.transportState = transportState;            
            state.G = transportState.G;
            
%             isDG = isfield(transportState, 'sdof');
            if model.isDG
                disc = model.transportModel.disc;
            end
            
            nPh = model.water + model.oil + model.gas;
            for vNo = 1:numel(varNames)
                vn  = varNames{vNo};
                dvn = dofVarNames{vNo};
                if isfield(transportState, vn)
                    if model.isDG
                        x = model.fineTransportModel.G.cells.centroids;
                        x = disc.transformCoords(x, p);
                        for phNo = 1:nPh
                            state.(vn)(:,phNo) = disc.evaluateSaturation(x, p, transportState.(dvn)(:, phNo), st);
                        end
                    else                        
                        state.(vn) = transportState.(vn)(p,:);
                    end
                end
            end
 
        end

    end
end

function [varNames, dofVarNames] = getTransportVarNames()
    
    varNames    = {'s'   };
    dofVarNames = {'sdof'};

end

function dgVarNames = getDGVarNames()

    dgVarNames = {'bfactor', 'degree', 'nDof'};

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

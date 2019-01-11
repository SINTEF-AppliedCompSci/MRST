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
        isReordering
        mappings
        computeCoarsePressure
        coarsePressureModel
    end
    
    methods
        function model = AdaptiveSequentialPressureTransportModel(pressureModel, transportModel, G0, varargin)
            model = model@SequentialPressureTransportModel(pressureModel, transportModel, varargin{:});
            model.G0 = G0;
            
            model.isReordering = isa(transportModel, 'ReorderingModel');
            tm = transportModel;
            if model.isReordering
                tm = transportModel.parent;
            end
            model.isDG = isa(tm, 'TransportBlackOilModelDG');
            
            model.fineTransportModel = tm;
            
            [model.mappings.newPartition, model.mappings.oldPartition] = deal(G0.partition);
            
            model.coarseTransportModel = model.upscaleTransportModelTPFA(G0.partition);
            if model.isReordering
                model.transportModel = transportModel;
                model.transportModel.parent = model.coarseTransportModel;
            else
                model.transportModel = model.coarseTransportModel;
            end
            
            model.coarseTransportModel.G.cells.refined = G0.cells.refined;
            model.coarseTransportModel.G.oldPartition = G0.partition;
            
            model.plotProgress = false;
            model.storeGrids   = true;
            
            model.computeCoarsePressure = false;
            model.coarsePressureModel = [];
            
            model = merge_options(model, varargin{:});
            
        end
        
        function [state, pressureReport, transportReport, pressure_ok, transport_ok, forceArg] = solvePressureTransport(model, state, state0, dt, drivingForces, iteration)
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
            pressureForceArg = forceArg;
            
            % First, solve the pressure using the pressure nonlinear
            % solver.
            [pressureState, pressureReport] = ...
                psolver.solveTimestep(state0, dt, model.pressureModel,...
                            'initialGuess', state, ...
                            forceArg{:});
            pressure_ok = pressureReport.Converged || psolver.continueOnFailure;
            
            if pressure_ok
                % If pressure converged, we proceed to solve the transport
                model.transportModel = state.transportModel;
                % Assemble transport residual
                transportState    = model.upscaleState(pressureState);
                transportState0   = model.upscaleState(state0);
                drivingForcesUpsc = model.mapForces(drivingForces);

                problem = model.transportModel.getEquations(transportState0, transportState, ...
                                   dt, drivingForcesUpsc, 'resOnly', true);
                iscell = strcmpi(problem.types, 'cell');
                residual = horzcat(problem.equations{iscell});
                if model.isReordering
                    G = model.transportModel.parent.G;
                else
                    G = model.transportModel.G;
                end
                if model.isDG
                    ix = rldecode((1:G.cells.num)', transportState.nDof, 1);
                    residual = accumarray(ix, abs(residual));
                    residual = residual./transportState.nDof;
                end

                % Make sure we dont refine already refined cells
                if isfield(model.transportModel.G.cells, 'refined') && 0
                    residual(model.transportModel.G.cells.refined) = 0;
                end
            
                tol = 1e-2;
                cells = abs(residual) > tol;
                tol   = 0.25*tol;
                cells(G.cells.refined & abs(residual) > tol) = true;
                [model, transportState, transportState0, drivingForces] ...
                    = model.refineTransportModel(cells, pressureState, state0, drivingForces);
                transportForces = drivingForces;
                
                if model.computeCoarsePressure
                    forceArg = model.pressureModel.getDrivingForces(transportForces);
                    [transportState, pr] = ...
                        psolver.solveTimestep(transportState0, dt, model.coarsePressureModel,...
                                    'initialGuess', transportState, ...
                                     forceArg{:});
                end

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
                
                [transportState, transportReport] = ...
                    tsolver.solveTimestep(transportState0, dt, model.transportModel,...
                                'initialGuess', transportState, ...
                                forceArg{:});
                transport_ok = transportReport.Converged;
                
                if model.plotProgress
                    if ishandle(1)
                        set(0, 'CurrentFigure', 1);
                    else
                        figure(1);
                    end
                    cla;
                    transportModel = model.transportModel;
                    if model.isReordering
                        transportModel = transportModel.parent;
                    end
                    G = transportModel.G;
                    plotCellData(G, transportState.s(:,1), 'edgecolor', 'none');
%                     plotGrid(G, 'facec', 'none', 'edgealpha', 0.3);
                    pink = [214, 154, 153]/255;
                    plotGrid(G, 'facec', 'none', 'edgecolor', pink);
                    
%                     cmap = mrstColormap('type', 'wateroil');
                    cmap = jet;
                    colormap(cmap)
                    axis equal tight
                    drawnow();
                end
                
                % Map transport quantities onto fine grid to prepare for
                % next pressure solve
                state = model.refineState(transportState, pressureState);
                state.transportModel = model.transportModel;
                forceArg = pressureForceArg;
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
            [transportModel, pressureModel] = model.upscaleTransportModelTPFA(mappings.newPartition);
            model.mappings = mappings;
            
            if isa(model.transportModel, 'ReorderingModel')
                model.transportModel.parent = transportModel;
            else
                model.transportModel = transportModel;
            end
            model.coarsePressureModel = pressureModel;
            transportState.pv = transportModel.operators.pv;
            transportState  = model.upscaleState(transportState);
            transportState0 = model.upscaleState(transportState0);
            
            if model.storeGrids
                transportState.G = transportModel.G;
            end

            forces = model.mapForces(forces);            
            
        end
        
        function forces = mapForces(model, forces)
            if model.isReordering
                G = model.transportModel.parent.G;
            else
                G = model.transportModel.G;
            end
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
            if isa(model.transportModel, 'ReorderingModel')
                transportModel = model.transportModel.parent;
                
            else
                transportModel = model.transportModel;
            end
            G  = transportModel.G;
            iscomp = isa(transportModel, 'ThreePhaseCompositionalModel');
            
            % Pore volumes
            pvf = model.fineTransportModel.operators.pv;
            pv  = transportModel.operators.pv;
            
            % Summation matrix
            S    = sparse(G.partition, (1:GF.cells.num)', 1);
%             pvbf          = bsxfun(@times, pvf, state.bfactor);
%             b = S*(pvf.*state.bfactor)./pv;
%             bf = state.bfactor;
            
            flowVarNames = getFlowVarNames();
            
            
            if iscomp
                % Compute coarse masses and mole fractions
                % 
                eos = transportModel.EOSModel;
                X = eos.getMassFraction(state.x);
                Y = eos.getMassFraction(state.y);
                assert(~model.water);
                mass = pvf.*(state.s(:, 1).*state.rho(:, 1).*X + state.s(:, 2).*state.rho(:, 2).*Y);
                
                mass_c = S*mass;
                
                state_c = struct();
                
                flowVarNames = {'pressure', 'T'};
                for fNo = 1:numel(flowVarNames)
                    vn = flowVarNames{fNo};
                    if isfield(state, vn)
                        state_c.(vn) = S*(pvf.*state.(vn))./(pv);
                    end
                end
                state_c.flux = state.flux;
                state_c.wellSol = state.wellSol;
                
                
                mass_fraction = mass_c./sum(mass_c, 2);
                state_c.components = eos.getMoleFraction(mass_fraction);
                state_c = transportModel.computeFlash(state_c, inf);
                
                getDensity = @(x, Z) eos.PropertyModel.computeDensity(state_c.pressure, x, Z, state_c.T, nan);
                
                rhoL = getDensity(state_c.x, state_c.Z_L);
                rhoV = getDensity(state_c.y, state_c.Z_V);
                
                % Total mass, integrated on fine-scale
                massT_f = sum(mass_c, 2);
                
                % Total mass, on new coarse scale
                massT_c = pv.*(rhoL.*state_c.s(:, 1) + rhoV.*state_c.s(:, 2));
                
                % Compensate for introduced mass-error
                sT = massT_f./massT_c;
                state_c.s = state_c.s.*sT;
                
                state = state_c;
                
%                 rmfields = {'x', 'y', 'mixing', 'switched', 'switchount', 'w', 'b', 'Z_L', 'Z_V'};
%                 for i = 1:numel(rmfields)
%                     if isfield(state, rmfields{i})
%                         state = rmfield(state, rmfields{i});
%                     end
%                 end
                
                
%                 disp('marker');
                
%                 mass = state.transportState.components.*state.transportState.rho.*sT;
%                 state.s = S*(pvPrev.*bPrev.*state.transportState.s)./(S*(pvPrev.*bPrev));

            else

                for fNo = 1:numel(flowVarNames)
                    vn = flowVarNames{fNo};
                    if isfield(state, vn)
                        state.(vn) = S*(pvf.*state.(vn))./(pv);
                    end
                end

                % Mapping from old to new coarse grid
                S = S*sparse((1:GF.cells.num)', state.transportState.G.partition, 1);

                pvPrev = state.transportState.pv;
    %             pvb = bsxfun(@times, pv, state.bfactor);

                bPrev  = state.transportState.bfactor;
                state.s = S*(pvPrev.*bPrev.*state.transportState.s)./(S*(pvPrev.*bPrev));

                ts = state.transportState;

                if model.isDG

                    disc = transportModel.disc;
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
                    state.degree = repmat(transportModel.disc.degree, G.cells.num, 1);            
    %             	state        = assignDofFromState(transportModel.disc, state);
                end
            end

            cFsign   = fineToCoarseSign(G);
            cFacesno = rldecode(1:G.faces.num, diff(G.faces.connPos), 2).';
            newFlux = zeros(G.faces.num, size(state.flux, 2));
            for i = 1:size(state.flux, 2)
                newFlux(:, i)   = accumarray(cFacesno, state.flux(G.faces.fconn, i) .* cFsign);
            end
            state.flux = newFlux;
        end
        
        function [transportModel, pressureModel] = upscaleTransportModelTPFA(model, partition, varargin)
            
            transportModel = upscaleModelTPFA(model.fineTransportModel, partition, varargin{:});
            pressureModel = [];
            if model.computeCoarsePressure
                pressureModel = upscaleModelTPFA(model.pressureModel, partition, varargin{:});
                pressureModel.extraStateOutput = true;
            end
            
            G = transportModel.G;
            if model.isDG
                G = coarsenCellDimensions(G);
            end
            G.cells.equal    = false;
            G.faces.equal    = false;
            G.cells.ghost    = false(G.cells.num, 1);
            G.cells.refined  = accumarray(G.partition, ones(G.parent.cells.num,1)) == 1;
            transportModel.G = G;
            
            if model.isDG
                
                if isCoarseGrid(model.transportModel.G) && 0
                    cub = model.updateCubature(transportModel);
                end
                
                disc = transportModel.disc;
                transportModel.disc = DGDiscretization(transportModel, ...
                      'degree'             , disc.degree             , ...
                      'basis'              , disc.basis.type         , ...
                      'useUnstructCubature', disc.useUnstructCubature, ...
                      'jumpTolerance'      , disc.jumpTolerance      , ...
                      'outTolerance'       , disc.outTolerance       , ...
                      'meanTolerance'      , disc.meanTolerance      , ...
                      'outLimiter'         , 'kill'                  , ...
                      'limitAfterConvergence', disc.limitAfterConvergence);
            end
            
        end
        
        function state = refineState(model, transportState, pressureState)
            
            [varNames, dofVarNames] = getTransportVarNames();
            if isa(model.transportModel, 'ReorderingModel')
                transportModel = model.transportModel.parent;
            else
                transportModel = model.transportModel;
            end
            iscomp = isa(transportModel, 'ThreePhaseCompositionalModel');
            p     = transportModel.G.partition;
            
            state = pressureState;
            state.wellSol = transportState.wellSol;
            
%             st.pvPrev = model.transportModel.operators.pv;
%             st.oldPartition = model.transportModel.G.partition;
            transportState.transportState = [];
                           
            state.transportState = transportState;
            if model.storeGrids
                state.G = transportState.G;
            end
            
            if iscomp
                % Downscale compositions
                eos = transportModel.EOSModel;
                state.components = transportState.components(p, :);
                state = transportModel.computeFlash(state, inf);
                
                
                getDensity = @(x, Z) eos.PropertyModel.computeDensity(state.pressure, x, Z, state.T, nan);
                rhoL = getDensity(state.x, state.Z_L);
                rhoV = getDensity(state.y, state.Z_V);
                pvf = model.pressureModel.operators.pv;
                pvc = model.transportModel.operators.pv;
                m = pvf.*(rhoL.*state.s(:, 1) + rhoV.*state.s(:, 2));
                
                % Integrated mass from the fine grid, on the coarse grid 
                massF = accumarray(p, m);
                % Coarse scale mass
                massC = pvc.*sum(transportState.rho.*transportState.s, 2);
                % These two should match, compensate by changing
                % saturations
                sT_c = massC./massF;
                state.s = state.s.*sT_c(p);
            else
                if model.isDG
                    disc = transportModel.disc;
                end
                nPh =   model.pressureModel.water ...
                      + model.pressureModel.oil   ...
                      + model.pressureModel.gas;
                for vNo = 1:numel(varNames)
                    vn  = varNames{vNo};
                    dvn = dofVarNames{vNo};
                    if isfield(transportState, vn)
                        if model.isDG
                            x = model.fineTransportModel.G.cells.centroids;
    %                         x = disc.transformCoords(x, p);
                            for phNo = 1:nPh
    %                             state.(vn)(:,phNo) = disc.evaluateSaturation(x, p, transportState.(dvn)(:, phNo), transportState);
                                state.(vn)(:,phNo) = disc.evaluateDGVariable(x, p, transportState, transportState.(dvn)(:, phNo));
                            end
                        else                        
                            state.(vn) = transportState.(vn)(p,:);
                        end
                    end
                end

                if model.isDG
                    dgVarNames = getDGVarNames();
                    for vNo = 1:numel(dgVarNames)
                        vn  = dgVarNames{vNo};
                        if isfield(transportState, vn)
                            state.(vn) = transportState.(vn)(p,:);
                        end
                    end
                end
            end
        end
        
        function cubature = updateCubature(model, transportModel)
           
            GF = model.fineTransportModel.G;
            Gprev = model.transportModel.G;
            Gnew  = transportModel.G;
            
            
            newp  = Gnew.partition;
            prevp = Gprev.partition;
            
%             coarseVolumeCub = model.coarseTransportModel.disc.volumeCubature;
            fineVolumeCub   = model.fineTransportModel.disc.volumeCubature;
            prevVolumeCub   = model.transportModel.disc.volumeCubature;
            
            
            numCellsNew  = accumarray(newp, ones(GF.cells.num,1));
            numCellsPrev = accumarray(prevp, ones(GF.cells.num,1));
            
            ix = numCellsNew(newp) == numCellsPrev(prevp);
            
            
%             new2fine = nan(max(newp),1);
%             new2fine(newp) = 1:GF.cells.num;
%             prev2fine = nan(max(prevp),1);
%             prev2fine(prevp) = 1:GF.cells.num;
%             
%             ix = (prev2fine == new2fine) & (numCellsNew == numCellsPrev);
%             mappings.newPartition - mappings.oldPartition;
%             
%             xPrev = prevVolumeCub.points;
            
            cubature = [];
            
            
            
            
        end

    end
end

function [varNames, dofVarNames] = getTransportVarNames()

%     flds = {'pressure', 's', 'x', 'y', 'components', 'K', 'b', 'mob', ...
%             'rho', 'flag', 'L', 'T', 'Z_L', 'Z_V', 'w', 'bfactor', ...
%             'w_p', 'dpressure', 'dpRel', 'dpAbs', 'switched', ...
%             'nDof', 'degree', 'sdof'};

    varNames    = {'s'   };
    dofVarNames = {'sdof'};

end

function varNames = getFlowVarNames()

%     flds = {'pressure', 's', 'x', 'y', 'components', 'K', 'b', 'mob', ...
%             'rho', 'flag', 'L', 'T', 'Z_L', 'Z_V', 'w', 'bfactor', ...
%             'w_p', 'dpressure', 'dpRel', 'dpAbs', 'switched', ...
%             'nDof', 'degree', 'sdof'};

    varNames    = {'pressure', 'mob', 'bfactor', 'dpRel', 'rho', 'flag', 'dpressure'};

end


function dgVarNames = getDGVarNames()

    dgVarNames = {'degree', 'nDof'};

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

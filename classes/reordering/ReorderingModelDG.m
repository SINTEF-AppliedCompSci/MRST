classdef ReorderingModelDG_ghost < ReorderingModel
    
    properties
        
        plotProgress
        plotAfterTimestep
        plotAfterCellSolve
        disc
        
    end
    
    methods
        function model = ReorderingModelDG_ghost(fullmodel, varargin)
            model = model@ReorderingModel(fullmodel);
            model.plotProgress       = false;
            model.plotAfterTimestep  = false;
            model.plotAfterCellSolve = false;
            model.disc               = model.parent.disc;
            model = merge_options(model, varargin{:});
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, forces, varargin)
           [problem, state] = model.parent.getEquations(state0, state, dt, forces, varargin{:});
        end
        
        function [state, report] = stepFunction(reorderModel, state, state0, dt,...
                                                drivingForces, linsolve, nonlinsolve,...
                                                iteration, varargin)
            % Prorp evaluation, operators etc. are handled by parent model.
            model = reorderModel.parent;
            % Get the ordering, account for any cycles due to gravity or
            % capillary forces
            [order, A] = getTopologicalFluxOrdering(model.G, state);
            % Check initial residual so we know when we stop. Note:
            % Assuming scaling and "simple" convergence here
            p0  = model.getEquations(state0, state, dt, drivingForces, ...
                                  'iteration', iteration, 'resOnly', true);
            bad = false(numel(p0.equations{1}), 1);
            for eqNo = 1:numel(p0.equations)
                sub      = abs(p0.equations{eqNo}) > model.nonlinearTolerance;
                bad(sub) = true;
            end
            % Map form bad dofs to bad cells
            ix  = rldecode((1:model.G.cells.num)', state.nDof, 1);
            bad = accumarray(ix, bad) > 0;
            ok  = ~bad;
            % If all cells are ok, make sure we at least perturb first one
            if all(ok)
                ok(order(1)) = false;
            end
            
            % Cells that has 
            changed        = ~ok;
            beenConsidered = false(model.G.cells.num, 1);
            
            % Do chunkSize cells at a time to avoid AD overhead slowdown
            chunkSize = reorderModel.chunkSize;
            
            % Book-keeping of cycles
            counts    = accumarray(order, 1);
            sumCounts = cumsum(counts);
            orderMax  = max(order);
            start     = 0;
            
            G            = model.G;
            state.solved = false(G.cells.num,1); % Logical vector of solved cells (for plotting)
            its          = zeros(G.cells.num,1); % Nonlinear iterations used for each cell
            while 1
                % Pick next group
                if start == 0
                    offset = 0;
                else
                    offset = sumCounts(start);
                end
                next = find(sumCounts >= offset + chunkSize, 1, 'first');
                if isempty(next)
                    next = orderMax;
                end
                cells = order > start & order <= next;
                
                % Assuming TPFA neighborship, we take the closest neighbors
                % to be changed
                touched = A'*changed;
                
                if any(cells) && any(touched(cells) | ~ok(cells))                    
                    % We have non-zero residual or changed values in our
                    % chunk, we need to do something
                    state.order  = order;
                    [state, rep] = reorderModel.performReordering(state, state0, dt, drivingForces, cells);
                    its(cells)   = rep.Iterations;
                    
                    % Mark cells as changed if any of the dofs have changed
                    % more than the nonlinear tolerance
                    chng = true(nnz(cells),model.disc.basis.nDof);
                    for dofNo = 1:model.disc.basis.nDof
                        ix  = model.disc.getDofIx(state , dofNo, cells, true);
                        ix0 = model.disc.getDofIx(state0, dofNo, cells, true);
                        ii  = ix > 0 & ix0 > 0;
                        dsdof = state.sdof(ix(ii),:) - state0.sdof(ix0(ii),:);
                        chng(ii, dofNo) = any(abs(dsdof) > model.nonlinearTolerance,2);
                    end
                    changed(cells) = any(chng,2);
                    
                    state0              = state;
                    state.solved(cells) = true;
                    
                    if reorderModel.plotProgress
                        
                        if ishandle(1)
                            set(0, 'CurrentFigure', 1);
                        else
                            fig          = figure(1);
                            fig.Position = [-1000, 0, 800, 800];
                            cmap         = mrstColormap();
                            colormap(cmap);
                        end
                        clf;
                        gr = [1,1,1]*0.3;
                        plotCellData(G, state.s(:,1), 'edgec', 'none');
                        plotGrid(G, 'facec', 'none', 'edgec', gr);
                        plotGrid(G, cells, 'facecolor', 'none', ...
                                           'edgecolor', 'r'   , ...
                                           'linewidth', 1.5   ); 
                        axis equal tight
                        drawnow();
                        
                        if 0
                        pth = fullfile(mrstPath('dg'), 'gifs', 'adaptive-reorder');
                        files = dir(pth);
                        names = {files.name};
                        names = names(3:end);
                        if isempty(names)
                            num = 1;
                        else
                            num = cellfun(@(n) str2double(n(5:end-4)), names);
                            num = max(num)+1;
                        end
                        
                        print(fullfile(pth, ['fig-', num2str(num)]), '-djpeg', '-r100');
                        end
                        
                    end
                    
                elseif all(ok(~touched))
                    % Everything downstream is already converged, quit
                    % early
                    break
                end
                % Set next strict lower bound to be current upper bound
                if next == orderMax
                    break
                end
                start = next;
                assert(~any(beenConsidered(cells)));
                beenConsidered(cells) = true;
            end
            report = model.makeStepReport('Converged'         , true, ...
                                          'ResidualsConverged', true, ...
                                          'Residuals'         , 0   );
            report.Iterations = its;

        end
        
        function [substate0, substate, submodel, subforces, mappings] = buildSubproblem(reordermodel, model, state, state0, forces, subs)
            % Build subproblem form full model
            
            % Fint active cells and faces for subgrid----------------------
            
            % Find active cells
            G  = model.G;
            Nc = G.cells.num;
            Nf = G.faces.num;
            submodel        = reordermodel.parent;
            keepCells       = false(Nc,1);
            keepCells(subs) = true;
            
            % Add all cells of wells intersecting the subset
            for wNo = 1:numel(forces.W)
                wc = forces.W(wNo).cells;
                if any(keepCells(wc))
                    keepCells(wc) = true;
                end
            end

            % Add ghost cells
            op0            = model.operators;
            active         = keepCells(op0.N);
            activeConn     = any(active,2);
            gc             = op0.N(activeConn,:);
            gc             = gc(~keepCells(gc));
            ghostCells     = false(Nc,1);
            ghostCells(gc) = true;
            keepCells(ghostCells) = true;

            GG = makeSubgrid(model.G, keepCells);
            
            % Active faces
            cells = find(keepCells);
            faces = G.cells.faces(mcolon(G.cells.facePos(cells)     , ...
                                         G.cells.facePos(cells+1)-1), 1);
            keepFaces        = false(G.faces.num,1);
            keepFaces(faces) = true;
            
            % Internal connections
            op0.internalConn(keepFaces,:);
            cells = find(keepCells & ~ghostCells);
            faces = G.cells.faces(mcolon(G.cells.facePos(cells)     , ...
                                         G.cells.facePos(cells+1)-1), 1);
            internalConn         = false(G.faces.num,1);
            internalConn(faces)  = true;
            bFaces               = boundaryFaces(G);
            internalConn(bFaces) = false;
            % Internal connections in subgrid
            internalConnParent = internalConn;
            internalConn       = internalConn(keepFaces);            
            
            % Mappings between sub and full grid
            c_o2n = zeros(Nc, 1);
            c_o2n(keepCells) = 1:nnz(keepCells);
            c_n2o = find(keepCells);
            cellMap = struct('keep'   , keepCells  , ...
                             'old2new', c_o2n      , ...
                             'new2old', c_n2o      , ...
                             'globalOrder'  , state.order, ...
                             'localOrder'   , c_o2n(state.order(keepCells(state.order))));
            
            f_o2n = zeros(Nf, 1);
            f_o2n(keepFaces) = 1:nnz(keepFaces);
            f_n2o = find(keepFaces);
            faceMap = struct('keep'   , keepFaces, ...
                             'old2new', f_o2n    , ...
                             'new2old', f_n2o    );
            
                         
            mappings = struct('cellMap', cellMap, ...
                              'faceMap', faceMap);
             
            nf = nnz(activeConn);
            nc = nnz(keepCells);
                          
            G.cells.ghost = false(G.cells.num,1);
            G.cells.ghost(ghostCells) = true;
            G.parent = G;
            
            G.cells.centroids = G.cells.centroids(keepCells, :);
            G.cells.volumes   = G.cells.volumes(keepCells);
%             G.cells.diameters = G.cells.diameters(keepCells);
            G.cells.dx        = G.cells.dx(keepCells,:);
            G.cells.num = nc;
            G.cells.ghost = false(G.cells.num,1);
            G.cells.ghost(mappings.cellMap.old2new(ghostCells)) = true;
            
            G.faces.centroids = G.faces.centroids(keepFaces, :);
            G.faces.areas = G.faces.areas(keepFaces);
            if isfield(G.faces, 'dx')
%                 G.faces.diameters = G.faces.diameters(keepFaces,:);
                G.faces.dx = G.faces.dx(keepFaces,:);
            end
            G.faces.num = nnz(keepFaces);
%             G.faces.neighbors = G.faces.neighbors(keepFaces,:);
%             G.faces.neighbors(G.faces.neighbors~=0) = cellMap.old2new(G.faces.neighbors(G.faces.neighbors~=0));

            G.mappings = mappings;
                          
            G.faces.neighbors = G.faces.neighbors(keepFaces,:);
            G.faces.neighbors(G.faces.neighbors~=0) = cellMap.old2new(G.faces.neighbors(G.faces.neighbors~=0));
            G.type = [G.type, {'subGrid'}];
%             internalConn = all(G.faces.neighbors > 0, 2);
%             internalConnParent = false(G.parent.faces.num,1);
%             internalConnParent(faceMap.new2old(internalConn)) = true;

            % Make operators for subgrid-----------------------------------
%             nf = nnz(activeConn);
%             nc = nnz(keepCells);
%             renumCell = zeros(1, nc);
%             renumCell(keepCells) = (1:nc);

            T     = op0.T(activeConn);
            T_all = op0.T_all(keepFaces);
            N = cellMap.old2new(op0.N(activeConn, :));
            M  = sparse((1:nf)'*[1 1], N, .5*ones(nf,2), nf, nc);
            C  = sparse( [(1:nf)'; (1:nf)'], N, ones(nf,1)*[1 -1], nf, nc);
            
            op = op0;
            op.C = C;
            op.M = M;
            upstr =  @(flag, x) faceUpstr(flag, x, N, [nf, nc]);
            op.faceUpstr = upstr;
            op.Div = @(x) C'*x;
            op.Grad = @(x) -C*x;
            op.faceAvg = @(x) M*x;
            op.splitFaceCellValue = @(operators, flag, x) splitFaceCellValue(operators, flag, x, [nf, nc]);
            
            disc = model.disc;
            vi = disc.velocityInterp;
            D = cell(1, G.griddim);
            for dNo = 1:G.griddim
                D{dNo} = vi.D{dNo}(keepCells, keepFaces);
            end
            vi.D = D;
            if G.griddim == 2
                vi.faceFlux2cellVelocity = @(v) [D{1}*v, D{2}*v];
            else
                vi.faceFlux2cellVelocity = @(v) [D{1}*v, D{2}*v, D{3}*v];
            end
            disc.velocityInterp = vi;
            
            op.T = T;
            op.T_all = T_all;
            op.N = N;
            op.pv = op0.pv(keepCells);
            op.internalConn = internalConn;
%             op.internalConn = true(nf, 1);
            
%             G.cells.ghost = false(G.cells.num,1);
%             G.cells.ghost(ghostCells) = true;
%             G.parent = G;
% 
%             G.cells.centroids = G.cells.centroids(keepCells, :);
%             G.cells.volumes   = G.cells.volumes(keepCells);
%             G.cells.diameters = G.cells.diameters(keepCells);
%             G.cells.dx        = G.cells.dx(keepCells,:);
%             G.cells.num = nc;
%             G.cells.ghost = false(G.cells.num,1);
%             G.cells.ghost(mappings.cellMap.old2new(ghostCells)) = true;
%             
%             G.faces.centroids = G.faces.centroids(keepFaces, :);
%             if isfield(G.faces, 'dx')
%                 G.faces.diameters = G.faces.diameters(keepFaces,:);
%                 G.faces.dx = G.faces.dx(keepFaces,:);
%             end
%             G.faces.num = nnz(keepFaces);
% %             G.faces.neighbors = G.faces.neighbors(keepFaces,:);
% %             G.faces.neighbors(G.faces.neighbors~=0) = cellMap.old2new(G.faces.neighbors(G.faces.neighbors~=0));
% 
%             G.mappings = mappings;

            % Copy state information---------------------------------------
            substate0 = state0;
            substate = state;
            flds = getCellFields();
            for fNo = 1:numel(flds)
                
                fn = flds{fNo};
                if isfield(substate, fn)
                    
                    ix = keepCells;
                    if size(state.(fn)(:,1),1) > Nc
                        ix = model.disc.getDofIx(state, Inf, cellMap.new2old);
                    end
                    substate.(fn) = state.(fn)(ix, :);
                    
                end

                if isfield(substate0, fn)
                    
                    ix0 = keepCells;
                    if size(state.(fn)(:,1),1) > Nc
                        ix0 = model.disc.getDofIx(state0, Inf, cellMap.new2old);
                    end
                    substate0.(fn) = state0.(fn)(ix0, :);
                end
                
            end
%             substate.flux = zeros(G.faces.num,2);
%             isbf = any(G.parent.faces.neighbors == 0,2);
%             substate.flux = state.flux(keepFaces & ~isbf,:);
%             substate.flux(
            substate.flux = state.flux(keepFaces,:);
%             isbf = any(G.faces.neighbors == 0,2);
%             substate.flux = substate.flux.*(internalConn | isbf);
%             substate.flux = state.flux(internalConnParent | isbf,:);

%             substate.dofPos  = substate.dofPos(:, keepCells);
%             substate0.dofPos = substate0.dofPos(:, keepCells);


            % Make well struct for submodel--------------------------------
            W = forces.W;

            keepWells = false(numel(W), 1);
            keep = keepCells & ~ghostCells;
            for wNo = 1:numel(W)
%                 keepw = keepCells(W(wNo).cells);
                keepw = keep(W(wNo).cells);
                if any(keepw)
                    assert(all(keepw));
                    W(wNo).cells = cellMap.old2new(W(wNo).cells);
                    keepWells(wNo) = true;
                end
            end
            W = W(keepWells);
            substate.wellSol  = substate.wellSol(keepWells);
            substate0.wellSol = substate0.wellSol(keepWells);
            
            % Make BC struct for submodel----------------------------------
            bc     = forces.bc;
            
            if ~isempty(bc)
                keepBC = any(bc.face == bFaces',2);
                bc.face  = faceMap.old2new(bc.face(keepBC));
                bc.type  = bc.type(keepBC);
                bc.value = bc.value(keepBC);
                bc.sat   = bc.sat(keepBC,:);
            end
            
            % Gather forces in force struct--------------------------------
            subforces    = forces;
            subforces.W  = W;
            subforces.bc = bc;
            
            % Make subrock
            rock = submodel.rock;
            rock.poro = rock.poro(keepCells);
            rock.perm = rock.perm(keepCells,:);
            submodel.rock = rock;
            
            % Make submodel------------------------------------------------
            submodel.operators = op;
            submodel.G = G;
            submodel.FacilityModel = FacilityModel(submodel);
            submodel.FacilityModel = submodel.FacilityModel.setupWells(subforces.W);
            
            disc.G = G;
            disc.internalConn = op.internalConn;
%             icp = false(G.parent.faces.num,1);
%             icp(ifaceList(active)) = op.internalConn;
            disc.internalConnParent = internalConnParent;
%             submodel.disc.internalConnParent = op0.internalConn;
            disc.N = N;
            submodel.disc = disc;
            
            substate = submodel.disc.updateDofPos(substate);
            substate0 = submodel.disc.updateDofPos(substate0);
            
            if 0
                figure(2)
                plotGrid(G.parent);
                plotGrid(G.parent, cellMap.keep, 'facec', 'r');
            end
            
%             dp = reshape((1:numel(bc.face)*disc.basis.nDof)', disc.basis.nDof, []);
%             nd = [ones(numel(faceBC),1); state.nDof(cellBCaux(:,1))];
%             bc.dofPos     = ones(
            
        end
         
        function [state, report] = performReordering(model, state, state0, dt, forces, subs)
            
            [substate0, substate, submodel, subforces] = model.buildSubproblem(model.parent, state, state0, forces, subs);
            forceArg = getDrivingForces(submodel, subforces);
            submodel.verbose = false;
            nls = NonLinearSolver();
            solveCell = find(submodel.G.mappings.cellMap.keep & ~submodel.G.parent.cells.ghost);
            if numel(solveCell) == 1
                fprintf('Solving cell %d ... ', solveCell);
            end
            [substate, report] = ...
                nls.solveTimestep(substate0, dt, submodel,...
                            'initialGuess', substate, ...
                            forceArg{:});
                        
            maps     = submodel.G.mappings.cellMap;
            keepFull = maps.keep & ~submodel.G.parent.cells.ghost;
            keepSub  = ~submodel.G.cells.ghost;
            
            state.s(keepFull, :)   = substate.s(keepSub,:);
            state.mob(keepFull, :) = substate.mob(keepSub,:);
            state.degree(keepFull) = substate.degree(keepSub);
            state.nDof(keepFull)   = substate.nDof(keepSub);
%             state.cfl(keepFull)    = substate.cfl(keepSub);
            
            nDof = model.parent.disc.basis.nDof;
            cells = (state.degree < model.parent.disc.degree) & keepFull;
            ix = model.parent.disc.getDofIx(state, 2:nDof, cells);
            state.sdof(ix,:) = [];    
% %                 ix = model.parent.disc.getDofIx(state, dofNo, keepFull);
%             ix = mcolon(state.nDof(keepFull)+1, model.parent.disc.basis.nDof);
%             state.sdof(ix,:) = [];
            
%             model.parent.disc.getDofIx(state, 
            
            state = model.parent.disc.updateDofPos(state);
            state = model.parent.disc.mapDofs(state, state0);
            
            ixFull = model.parent.disc.getDofIx(state, Inf, keepFull);
            ixSub  = submodel.disc.getDofIx(substate, Inf, keepSub);
            state.sdof(ixFull,:) = substate.sdof(ixSub,:);
            
            if model.plotAfterCellSolve
                        
                if ishandle(1)
                    set(0, 'CurrentFigure', 1);
                else
                    figure(1);
                end
                clf;
                nr = 1; nc = 2; G = model.G;

                subplot(nr,nc,1);
                plotCellData(G, state.s(:,1), 'edgec', 'none');
                colorbar();
                axis equal tight
                
                subplot(nr,nc,2);
                plotGrid(G, 'facec', 'none', 'edgec', [1,1,1]*0.6);
                plotGrid(G, keepFull, 'edgec', [1,1,1]*0.6);
                axis equal tight

                drawnow();

            end
            
        end
        
        function bc = setBoundaryConditions(model, bc, option, isType, pressure, volFluxRes)
            
            assert
            
            switch lower(option)
                case 'volume_flux'
                    type = 'rflux';
                    vf = volFluxRes(isType, :);
                    v = sum(vf, 2);
                    s = bsxfun(@rdivide, abs(vf), sum(abs(vf), 2));
%                 case 'mass_flux'
%                     type = 'flux';
%                     vf = volFluxSurf(isType, :);
%                     v = sum(vf, 2);
%                     s = bsxfun(@rdivide, abs(vf), sum(abs(vf), 2));
                case 'pressure'
                    type = 'pressure';
                    v = pressure(isType);
                    vf = volFluxRes(isType, :);
                    s = bsxfun(@rdivide, abs(vf), sum(abs(vf), 2));
                case 'component_flux'
                    disp('hello_world');
                otherwise
                    error('Unknown treatment');
            end
            [bc.type{isType}] = deal(type);
            bc.value(isType, :) = v;
            bc.sat(isType, :) = s;
        end
        
    end
end

function flds = getCellFields()
%     flds = {'pressure', 's', 'x', 'y', 'components', 'K', 'b', 'mob', ...
%             'rho', 'flag', 'L', 'T', 'Z_L', 'Z_V', 'w', 'bfactor', ...
%             'w_p', 'dpressure', 'dpRel', 'dpAbs', 'switched', ...
%             'nDof', 'degree', 'sdof'};
    flds = {'pressure', 's', 'x', 'y', 'components', 'K', 'b', 'mob', ...
            'rho', 'flag', 'L', 'T', 'Z_L', 'Z_V', 'w', 'bfactor', ...
            'w_p', 'dpressure', 'dpRel', 'dpAbs', 'switched', ...
            'nDof', 'degree', 'sdof'};
end
    
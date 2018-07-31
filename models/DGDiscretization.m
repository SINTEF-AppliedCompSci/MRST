classdef DGDiscretization < HyperbolicDiscretization
    
    properties

        dim
        degree
        basis
        volumeCubature
        surfaceCubature
        useUnstructCubature
        jumpTolerance
        outTolerance
        meanTolerance
        Gfull
        
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function disc = DGDiscretization(model, dim, varargin)
            
            disc = disc@HyperbolicDiscretization(model);
            
            disc.dim    = dim;
            disc.degree = 1;
            disc.basis  = 'legendre';
            disc.jumpTolerance = 0.2;
            disc.outTolerance = 1e-4;
            disc.meanTolerance = 1e-4;
            disc.useUnstructCubature = false;
            
            disc        = merge_options(disc, varargin{:});
            
            disc.basis  = dgBasis(dim, disc.degree, disc.basis);
            
            G = disc.G;
            disc.Gfull = G;
            if G.griddim == 2
                if disc.degree == 0 || disc.useUnstructCubature
                    disc.volumeCubature = Unstruct2DCubature(disc.G, disc.degree + 1, disc.internalConn);
                else
                    disc.volumeCubature  = TriangleCubature(disc.G, disc.degree + 1, disc.internalConn);
                end
                disc.surfaceCubature  = LineCubature(disc.G, disc.degree + 1, disc.internalConn);
            else
                if disc.degree == 0 || disc.useUnstructCubature
                    disc.volumeCubature = Unstruct3DCubature(disc.G, disc.degree + 1, disc.internalConn);
                    disc.surfaceCubature = Unstruct2DCubature(disc.G, disc.degree + 1, disc.internalConn);
%                     disc.surfaceCubature = TriangleCubature(disc.G, disc.degree + 1, disc.internalConn);
                else
                    disc.volumeCubature  = TetrahedronCubature(disc.G, disc.degree + 1, disc.internalConn);
                    disc.surfaceCubature = TriangleCubature(disc.G, disc.degree + 1, disc.internalConn);
                end
            end
%             [W, x, w, ii, jj, cellNo] = disc.volumeCubature.getCubature((1:G.cells.num)', 'cell');
%             disc.volumeCubature.W = W;
%             [W, x, w, ii, jj, cellNo] = disc.surfaceCubature.getCubature((1:G.cells.num)', 'cellsurface');
%             disc.surfaceCubature.W = W;
            
        end
        
        %-----------------------------------------------------------------%        
        function state = assignDofFromState(disc, state)

            state.degree = repmat(disc.degree, disc.G.cells.num, 1);
            state = disc.updateDisc(state);
            
            state.nDof = disc.getnDof(state);
            sdof = zeros(sum(state.nDof), size(state.s,2));
            
            ix = disc.getDofIx(state, 1);
            sdof(ix, :) = state.s;

            state.sdof = sdof;

        end
        
        %-----------------------------------------------------------------%
        function state = updateDisc(disc, state)
            
%             Nc = disc.G.cells.num;
%             if isfield(state, 'mappings')
%                 Nc = nnz(state.mappings.cellMap.keepCells);
%             end
%             dp = reshape((1:disc.G.cells.num*disc.basis.nDof)', disc.basis.nDof, []);
            dp = reshape((1:numel(state.cells)*disc.basis.nDof)', disc.basis.nDof, []);

%             dp = reshape((1:Nc*disc.basis.nDof)', disc.basis.nDof, []);
            
            if nargin == 1
                nd = repmat(disc.basis.nDof, numel(state.cells), 1);
            else
                nd = disc.getnDof(state);
                subt = cumsum([0; disc.basis.nDof - nd(1:end-1)]);
                [ii, jj, v] = find(dp);
                
                if size(ii,1) == 1, ii = ii'; end
                if size(jj,1) == 1, jj = jj'; end
                if size(v ,1) == 1, v  =  v'; end

                cnDof = cumsum(nd);

                v = v - subt(jj);
                v(v > cnDof(jj)) = 0;
                dp = full(sparse(ii, jj, v));
            end
            
            state.nDof   = nd;
            state.dofPos = dp;
            
%             disc.state = state;
            
        end
        
        %-----------------------------------------------------------------%
        function [xhat, translation, scaling] = transformCoords(disc, x, cells, inverse)
            
            G = disc.G;
            translation = -G.cells.centroids(cells,:);
            if isfield(G.cells, 'dx')
                scaling = 1./(G.cells.dx(cells,:)/2);
            else
                scaling = 1./(G.cells.diameters(cells)/(2*sqrt(G.griddim)));
            end
            
            if nargin < 4 || ~inverse
                xhat = (x + translation).*scaling;
                xhat = xhat(:, 1:disc.dim);
                scaling = scaling(:, 1:disc.dim);
                translation = translation(:, 1:disc.dim);
            else
                xhat = x./scaling - translation;
            end
               
        end
        
        %-----------------------------------------------------------------%
        function v = trimValues(disc, v)
            
            tol = eps(mean(disc.G.cells.volumes));
            tol = -inf;
%             tol = 1e-7;
            ix = abs(v) < tol;
            if isa(v, 'ADI')
                v.val(ix) = 0;
            else
                v(ix) = 0;
            end
            
        end
        
        %-----------------------------------------------------------------%
        function ix = getDofIx(disc, state, dofNo, cells, includezero)

            G = disc.G;
            if nargin == 2
                cells = 1:disc.G.cells.num;
                dofNo = 1:disc.basis.nDof;
            elseif nargin == 3
                cells = 1:G.cells.num;
            else
                if isempty(dofNo)
                    dofNo = 1:disc.basis.nDof;
                end
            end
            
            ix = state.dofPos(dofNo, cells);
            ix = ix(:);
            
            if nargin < 5
                includezero = false;
            end
            if ~includezero
                ix(ix == 0) = [];
            end
              
        end
        
        %-----------------------------------------------------------------%
        function nDof = getnDof(disc, state)
            
            if disc.degree < 0
                nDof = 0;
            else
                nDof = factorial(state.degree + disc.dim)...
                       ./(factorial(disc.dim).*factorial(state.degree));
            end
            
        end
        
        %-----------------------------------------------------------------%
        function state = mapDofs(disc, state, state_prev)
            
            state = disc.updateDisc(state);
            
            if all(state.nDof == state_prev.nDof)
                return
            else
                sdof = zeros(sum(state.nDof), size(state_prev.sdof,2));
                for dofNo = 1:disc.basis.nDof
                    ix      = disc.getDofIx(state, dofNo, (1:disc.G.cells.num)', true);
                    ix_prev = disc.getDofIx(state_prev, dofNo, (1:disc.G.cells.num)', true);
                    
                    sdof(ix(ix_prev > 0),:) = state_prev.sdof(ix_prev(ix_prev > 0),:);
                    
                end
                state.sdof = sdof;
            end

        end
        
        %-----------------------------------------------------------------%
        function s = evaluateSaturation(disc, x, cells, dof, state)
            
            psi     = disc.basis.psi;
            nDof    = state.nDof;
            nDofMax = disc.basis.nDof;
            
            ix = disc.getDofIx(state, 1, cells);
            s = dof(ix).*0;
            for dofNo = 1:nDofMax
                keep = nDof(cells) >= dofNo;
                ix = disc.getDofIx(state, dofNo, cells(keep));              
                s(keep) = s(keep) + dof(ix).*psi{dofNo}(x(keep,:));
            end
            
        end
        
        %-----------------------------------------------------------------%
        function state = getCellSaturation(disc, state)
            
            
%             [W, x, cellNo, faceNo] = disc.getCubature((1:disc.G.cells.num)', 'volume');
            [W, x, cellNo, faceNo] = disc.getCubature(state.cells, 'volume');
%             getCubature(disc, cells, type)
%             [x, w, nq, ii, jj, cellNo] = makeCellIntegrator(disc.G, (1:disc.G.cells.num)', max(disc.degree), 'volume');
%             W = sparse(ii, jj, w);

            x = disc.transformCoords(x, cellNo);
            
            sdof = state.sdof;
            nPh = size(sdof,2);
            s = zeros(disc.G.cells.num, nPh);
            for phNo = 1:nPh
                s(:,phNo) = (W*disc.evaluateSaturation(x, cellNo, sdof(:,phNo), state))./disc.G.cells.volumes(state.cells);
            end
            
            state.s = s;
            
        end
        
        %-----------------------------------------------------------------%
        function I = cellInt(disc, integrand, f, cells, sdof, sdof0, state, state0)
        
            G    = disc.G;
            psi  = disc.basis.psi;
            grad_psi = disc.basis.grad_psi;
            nDof = state.nDof;
            nDofMax = disc.basis.nDof;
            
            [W, x, cellNo, faceNo] = disc.getCubature(cells, 'volume');
            [x, ~, scaling] = disc.transformCoords(x, cellNo);

            s  = disc.evaluateSaturation(x, cellNo , sdof , state );
            s0 = disc.evaluateSaturation(x, cellNo, sdof0, state0);
            f = f(s, 1-s, cellNo, cellNo);
            
            I = integrand(sdof, sdof, sdof, ones(numel(double(sdof)), 1), 1, ones(1, disc.dim)).*0;
            
            for dofNo = 1:nDofMax
                
                keepCells = nDof(cells) >= dofNo;
                
                if any(keepCells)
                
                    ix = disc.getDofIx(state, dofNo, cells(keepCells));
                    i  = W*integrand(s, s0, f, cellNo, psi{dofNo}(x), grad_psi{dofNo}(x).*scaling);
                    I(ix) = i(keepCells);
                    
                elseif numel(cells) == disc.G.cells.num
                    
                    warning('No cells with %d dofs', dofNo);
                    
                end
                
            end
            
            I = disc.trimValues(I);
            
        end
        
        %-----------------------------------------------------------------%
        function I = faceFluxInt(disc, integrand, f, cells, sdof, state, T, vT, G, mob)
            
            g       = disc.G;
            psi     = disc.basis.psi;
            nDof    = state.nDof;
            nDofMax = disc.basis.nDof;
            
            [W, x, cellNo, faceNo] = disc.getCubature(state.cells(cells), 'surface');
            nPh = numel(G);
            
            [flag_v, flag_G, cL, cR] = disc.getUpstreamCell(faceNo, x, T, vT, G, mob, sdof, state);
            
            [upCells_v, upCells_G] = deal(repmat({cR}, 1, nPh));
            [x_v, s_v, x_G, s_G] = deal(cell(1, nPh));
            for phNo = 1:nPh
                
                upCells_v{phNo}(flag_v(:,phNo)) = cL(flag_v(:,phNo));
                x_v{phNo} = disc.transformCoords(x, upCells_v{phNo});
                s_v{phNo} = disc.evaluateSaturation(x_v{phNo}, upCells_v{phNo}, sdof, state);
                
                upCells_G{phNo}(flag_G(:,phNo)) = cL(flag_G(:,phNo));
                x_G{phNo} = disc.transformCoords(x, upCells_G{phNo});
                s_G{phNo} = disc.evaluateSaturation(x_G{phNo}, upCells_G{phNo}, sdof, state);
                
            end
            
            f_v = f(s_v{1}, 1-s_v{2}, upCells_v{1}, upCells_v{2});
            f_G = f(s_G{1}, 1-s_G{2}, upCells_G{1}, upCells_G{2});
            
            [x_c, ~, ~] = disc.transformCoords(x, cellNo);

            I = integrand(sdof, sdof, sdof, sdof, 1, 1, 1, 1, 1).*0;
            
            for dofNo = 1:nDofMax
                
                keepCells = nDof(cells) >= dofNo;
                
                if any(keepCells)
                    
                    ix = disc.getDofIx(state, dofNo, cells(keepCells)');
                    i  = W*integrand(f_v, f_G, s_v{1}, s_G{1}, cellNo, upCells_v{1}, upCells_G{1}, faceNo, psi{dofNo}(x_c));
                    I(ix) = i(keepCells);
                    
                elseif numel(cells) == disc.G.cells.num
                    
                    warning('No cells with %d dofs', dofNo);
                    
                end
                
            end
            
            I = disc.trimValues(I);
            
        end
        
        %-----------------------------------------------------------------%
        function I = faceFluxIntBC(disc, integrand, f, sdof, state, bc)
            
            G       = disc.G;
            psi     = disc.basis.psi;
            nDof    = state.nDof;
            nDofMax = disc.basis.nDof;
            
            isFlux = strcmpi(bc.type, 'flux');
            faces = bc.face(isFlux);
            
            cells = bc.cell(:,1);
            
            [W, x, cellNo, faceNo] = disc.getCubature(faces, 'face');
%             
%             [W, x, cellNo, faceNo] = disc.getCubature(state.cells(cells), 'surface');
%             nPh = numel(G);
%             
%             [flag_v, flag_G, cL, cR] = disc.getUpstreamCell(faceNo, x, T, vT, G, mob, sdof, state);
%             
%             [upCells_v, upCells_G] = deal(repmat({cR}, 1, nPh));
%             [x_v, s_v, x_G, s_G] = deal(cell(1, nPh));
%             for phNo = 1:nPh
%                 
%                 upCells_v{phNo}(flag_v(:,phNo)) = cL(flag_v(:,phNo));
%                 x_v{phNo} = disc.transformCoords(x, upCells_v{phNo});
%                 s_v{phNo} = disc.evaluateSaturation(x_v{phNo}, upCells_v{phNo}, sdof, state);
%                 
%                 upCells_G{phNo}(flag_G(:,phNo)) = cL(flag_G(:,phNo));
%                 x_G{phNo} = disc.transformCoords(x, upCells_G{phNo});
%                 s_G{phNo} = disc.evaluateSaturation(x_G{phNo}, upCells_G{phNo}, sdof, state);
%                 
%             end
%             
%             f_v = f(s_v{1}, 1-s_v{2}, upCells_v{1}, upCells_v{2});
%             f_G = f(s_G{1}, 1-s_G{2}, upCells_G{1}, upCells_G{2});
%             
%             [x_c, ~, ~] = disc.transformCoords(x, cellNo);
            
            
%             [flag_v, flag_G, cL, cR] = disc.getUpstreamCell(
            

            cellNo = sum(G.faces.neighbors(faceNo,:),2);
            
            sgn = 1 - 2*(G.faces.neighbors(faces, 1) ~= cells);
            W = W.*sgn;
            numFlux = nnz(isFlux);
            
            [xR, ~, ~] = disc.transformCoords(x, cellNo);
            sR = disc.evaluateSaturation(xR, cellNo, sdof, state);
            
            all2BC = nan(G.faces.num,1);
            all2BC(faces) = 1:numFlux;
            locFaceNo = all2BC(faceNo);
            
            sL = bc.sat(isFlux,1);
            sL = sL(locFaceNo);
            
            isInj = bc.value(isFlux) > 0;
            s = sL.*isInj(locFaceNo) + sR.*(~isInj(locFaceNo));
            
            
            f = f(s, 1-s, cellNo, cellNo);
            
            I = integrand(sdof, sdof, 1, 1, 1).*0;
            
            for dofNo = 1:nDofMax
                
                keepCells = nDof(cells) >= dofNo;
                
                if any(keepCells)
                    
                    ix = disc.getDofIx(state, dofNo, cells(keepCells)');
                    i  = W*integrand(f, s, cellNo, faceNo, psi{dofNo}(xR));
                    I(ix) = i(keepCells);
                    
                end
                
            end
            
            I = disc.trimValues(I);
            
        end
        
        function [flag_v, flag_G, cL, cR] = getUpstreamCell(disc, faces, x, T, vT, G, mob, sdof, state)
            
            g = disc.G;
            all2int = zeros(g.faces.num,1);
            all2int(disc.internalConn) = 1:nnz(disc.internalConn);
            ix = all2int(faces);
            
            cL = disc.N(ix,1);
            cR = disc.N(ix,2);
            T = T(ix);
            vT = vT(faces);
%             G = cellfun(@(g) g(ix), G, 'unif', false);
            G = cellfun(@(g) g(faces), G, 'unif', false);
            
            [xL, ~, ~] = disc.transformCoords(x, cL);
            [xR, ~, ~] = disc.transformCoords(x, cR);
            
            sL = disc.evaluateSaturation(xL, cL, sdof, state);
            sR = disc.evaluateSaturation(xR, cR, sdof, state);
            mob{1} = mob{1}([sL; sR], [cL; cR]);
            mob{2} = mob{2}(1-[sL; sR], [cL; cR]);
            
            N = [1:numel(ix); numel(ix)+1:2*numel(ix)]';
            
            upw = @(flag, x)faceUpstr(flag, x, N, [size(N,1), max(max(N))]);
            [flag_v, flag_G] = getSaturationUpwind('potential', 0, G, vT, T, mob, upw);
            
        end
        
        %-----------------------------------------------------------------%
        function [W, x, cellNo, faceNo] = getCubature(disc, elements, type)
            
            if size(elements,1) == 1, elements = elements'; end
            
            switch type
                case 'volume'
                    
                    cubature = disc.volumeCubature;
                    ix = mcolon(cubature.parentPos(elements), cubature.parentPos(elements+1)-1);
%                     ixf = ix;
                    nq = diff(cubature.parentPos);
                    nq = nq(elements);
                    
                    cellNo = rldecode(elements, nq, 1);
                    faceNo = [];

                    sgn   = ones(numel(ix),1);
                    
                case 'face'
                    
                    cubature = disc.surfaceCubature;
                    nqf = diff(cubature.parentPos);
                    faceNo = rldecode(elements, nqf(elements), 1);
                    ix = mcolon(cubature.parentPos(elements), cubature.parentPos(elements+1)-1);
                    nq = nqf(elements);
                    cellNo = [];
                    sgn = 1;
                    
                case 'surface'
                    
                    cubature = disc.surfaceCubature;
                    G = disc.G;
                    faces = G.cells.faces(mcolon(G.cells.facePos(elements), G.cells.facePos(elements+1)-1));
                    ncf = accumarray(rldecode((1:G.cells.num)', diff(G.cells.facePos), 1), disc.internalConn(G.cells.faces(:,1)));
                    faces = faces(disc.internalConn(faces));
                    if size(faces,1) == 1, faces = faces'; end
                    ix = mcolon(cubature.parentPos(faces), cubature.parentPos(faces+1)-1);
                    
                    nqf = diff(cubature.parentPos);
                    nq = accumarray(rldecode(elements, ncf(elements), 1), nqf(faces));
                    
                    cellNo = rldecode(elements, nq);
                    faceNo = rldecode(faces, nqf(faces), 1);
                    sgn   = 1 - 2*(G.faces.neighbors(faceNo,1) ~= cellNo);
                    
%                     ixf = nan(numel(ix),1);
%                     ixf(disc.internalConn(faceNo)) = 1:nnz(disc.internalConn(faceNo));
                    
            end
            
            x = cubature.points(ix, :);
            w = cubature.weights(ix).*sgn;
            [ii, jj] = blockDiagIndex(ones(numel(elements),1), nq);
            W = sparse(ii, jj, w);
%             W = cubature.W(cells, ixf);
            
        end
        
        %-----------------------------------------------------------------%
        function [smin, smax] = getMinMaxSaturation(disc, state)
            
            G = disc.G;
            
            faces = G.cells.faces(:,1);
            nodes = G.faces.nodes(mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces+1)-1));
            
            nfn = diff(G.faces.nodePos);
            ncn = accumarray(rldecode((1:G.cells.num)', diff(G.cells.facePos), 1), nfn(faces));
            
            x = G.nodes.coords(nodes,:);
            
            cells = rldecode((1:G.cells.num)', ncn, 1);
            x = disc.transformCoords(x, cells);
            
            s = disc.evaluateSaturation(x, cells, state.sdof(:,1), state);
            
            jj = rldecode((1:G.cells.num)', ncn, 1);
            s = sparse(jj, (1:numel(s))', s);
            
            smax = full(max(s, [], 2));
            smin = full(min(s, [], 2));
            
        end
        
        %-----------------------------------------------------------------%
        function jump = getInterfaceJumps(disc, state)
            
            G = disc.G;
            
            faces = G.cells.faces(:,1);
            isbf  = any(G.faces.neighbors(faces,:) == 0,2);
            faces = faces(~isbf);
            
            s     = @(x, c) disc.evaluateSaturation(x, c, state.sdof(:,1), state);
            
            cells = rldecode((1:G.cells.num)', diff(G.cells.facePos), 1);
            nbf   = accumarray(cells, isbf);
            cells = rldecode((1:G.cells.num)', diff(G.cells.facePos) - nbf, 1);
            
            xf = G.faces.centroids(faces,:);
            
            c_l  = G.faces.neighbors(faces,1);
            xf_l = disc.transformCoords(xf, c_l);
            c_r  = G.faces.neighbors(faces,2);
            xf_r = disc.transformCoords(xf, c_r);
            
            
            jump = abs(s(xf_l, c_l) - s(xf_r, c_r));
            
            jump = accumarray(cells, jump > disc.jumpTolerance) > 0;
            
            
        end
        
        %-----------------------------------------------------------------%
        function state = limiter(disc, state)
            
            G = disc.G;
            nDofMax = disc.basis.nDof;

            [smin, smax] = disc.getMinMaxSaturation(state);
            
            jump = disc.getInterfaceJumps(state);
            
            over  = smax > 1 + disc.outTolerance;
            under = smin < 0 - disc.outTolerance;
            outside = over | under;
            
            s = disc.getCellSaturation(state);
            s = s.s(:,1);
            mover  = s > 1 + disc.meanTolerance;
            munder = s < 0 - disc.meanTolerance;
            moutside = mover | munder;
            
            bad     = outside | moutside | jump;

            state.outside = outside;
            state.moutside = moutside;
            state.jump = jump;
            
%             state0 = state;
%             state.degree = repmat(disc.degree, G.cells.num, 1);
%             state = disc.mapDofs(state, state0);
%             state = disc.updateDisc(state);
            
            sdof = state.sdof;
           
            if any(bad)
                
                cells = find(bad);
                
                state.degree(cells) = 0;
                
%                 disc.dofPos = disc.updateDofPos();

                ix = disc.getDofIx(state, 1, cells);
                sdof(ix,:) = min(max(sdof(ix,:), 0), 1);
                sdof(ix,:) = sdof(ix,:)./sum(sdof(ix,:),2);
                
                state.s(cells,:) = sdof(ix,:);
                
                ix = disc.getDofIx(state, 2:nDofMax, cells);
                sdof(ix,:) = [];
                
                state.sdof = sdof;
                
                
%                 sdof = sdof(:,1);
%                 cells_o = find(over);
%                 cells_u = find(under);
%                 
%                 % Reduce to first-order
%                 ix = disc.getDofIx(state, (1 + G.griddim + 1):nDofMax, cells);
%                 sdof(ix) = 0;
%                 
%                 % Reduce linear dofs so that saturaion is within [0,1]
% %                 ix0 = disc.getDofIx(1                , cells);
%                 for dNo = 1:G.griddim
%                     
%                     ix0 = disc.getDofIx(state, 1, cells_u);
%                     ix  = disc.getDofIx(state, dNo+1, cells_u);
%                     sdof(ix) = sign(sdof(ix)).*min(abs(sdof(ix0)), abs(sdof(ix)));
%                     
%                     ix0 = disc.getDofIx(state, 1, cells_o);
%                     ix  = disc.getDofIx(state, dNo+1, cells_o);
%                     sdof(ix) = sign(sdof(ix)).*(1 - sdof(ix0));
%                     
%                 end
%                 
%                 state.sdof(:,1) = sdof;
%                 state.sdof(:,2) = -sdof;
%                 
%                 ix0 = disc.getDofIx(state, 1, (1:G.cells.num)');
%                 state.sdof(ix0,2) = 1 - state.sdof(ix0,1);
%                 
%                 [smin, smax] = disc.getMinMaxSaturation(state);
%                 
%                 s = disc.getCellSaturation(state);
                
            end
            
        end
        
        %-----------------------------------------------------------------%
        function v = faceFlux2cellVelocity(disc, model, faceFlux)
            
            
            
        end
        
        %-----------------------------------------------------------------%
        function plotCellSaturation(disc, state, cellNo)
            
            G = disc.G;
            
            faces = G.cells.faces(G.cells.facePos(cellNo):G.cells.facePos(cellNo+1)-1);
            nodes = G.faces.nodes(mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces+1)-1));
            nodes = reshape(nodes, 2, [])';
            
            swap = G.faces.neighbors(faces,1) ~= cellNo;

            nodes(swap,:) = nodes(swap, [2,1]); nodes = nodes(:,1);
            
            x = G.nodes.coords(nodes,:);
            
            x = disc.transformCoords(x, cellNo);
            
            
            if disc.dim == 1
                
                n = 100;
                xx = linspace(-1,1,n)';
                xk = xx;
                
            elseif disc.dim == 2
                
                n = 10; 
                xx = linspace(-1, 1, n);
                [xx, yy] = ndgrid(xx);
                xx = [xx(:), yy(:)];
                
                [in, on] = inpolygon(xx(:,1), xx(:,2), x(:,1), x(:,2));
                keep = in;
                xk = xx(keep,:);
                
            elseif disc.dim == 3
                
            end
                
            cellNo = repmat(cellNo, size(xk,1), 1);
            s = disc.evaluateSaturation(xk, cellNo, state.sdof, state);
            
            if disc.dim > 1
                s = scatteredInterpolant(xk, s);
                s = reshape(s(xx(:,1), xx(:,2)), n,n);
                surf(s);
            else
                plot(xx, s);
            end
            
        end
        
    end
        
end
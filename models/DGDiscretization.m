classdef DGDiscretization < HyperbolicDiscretization% < WENODiscretization
    
    properties
%         G
        dim
        degree
        basis
%         dofPos
%         nDof
%         state
%         limiter
%         cellIntegrator
%         faceIntegrator
    end
    
    methods
        %-----------------------------------------------------------------%
        function disc = DGDiscretization(model, dim, varargin)
            
%             disc = disc@WENODiscretization(model, dim, 'interpolateReference', false);
            disc = disc@HyperbolicDiscretization(model);
            
%             disc.G = model.G;
            disc.dim = dim;
            
            disc.degree = 1;
            disc.basis  = 'legendre';
            disc        = merge_options(disc, varargin{:});
            
            disc.basis  = dgBasis(dim, disc.degree, disc.basis);
%             disc.degree = repmat(disc.degree, disc.G.cells.num, 1);
%             [~, disc]        = disc.updateDisc();
            
%             disc.dofPos = reshape((1:disc.G.cells.num*disc.basis.nDof)', disc.basis.nDof, []);
            
        end
        
        %-----------------------------------------------------------------%        
        function state = assignDofFromState(disc, state)

            state.degree = repmat(disc.degree, disc.G.cells.num, 1);
            
%             if ~isfield(state, 'dofPos')
                state = disc.updateDisc(state);
%             end
            
            state.nDof = disc.getnDof(state);
            sdof = zeros(sum(state.nDof), size(state.s,2));
            
            ix = disc.getDofIx(state, 1);
            sdof(ix, :) = state.s;

            state.sdof = sdof;

        end
        
        %-----------------------------------------------------------------%
        function state = updateDisc(disc, state)
            
            dp = reshape((1:disc.G.cells.num*disc.basis.nDof)', disc.basis.nDof, []);
            
            if nargin == 1
                nd = repmat(disc.basis.nDof, disc.G.cells.num, 1);
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
        function [xhat, translation, scaling] = transformCoords(disc, x, cells)
            
            G = disc.G;
            translation = -G.cells.centroids(cells,:);
            if isfield(G.cells, 'dx')
                scaling = 1./(G.cells.dx(cells,:)/2);
            else
                scaling = 1./(G.cells.diameters(cells)/(2*sqrt(G.griddim)));
            end
            
            xhat = (x + translation).*scaling;
            
            xhat = xhat(:, 1:disc.dim);
            scaling = scaling(:, 1:disc.dim);
            translation = translation(:, 1:disc.dim);
               
        end
        
        %-----------------------------------------------------------------%
        function v = trimValues(disc, v)
            
            tol = eps(mean(disc.G.cells.volumes));
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
%             else
%                 ix(ix == 0) = 1;
%             end
              
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
%             s = double2ADI(zeros(numel(cells), 1), dof(ix));
            s = dof(ix).*0;
            for dofNo = 1:nDofMax
                keep = nDof(cells) >= dofNo;
                ix = disc.getDofIx(state, dofNo, cells(keep));
                
%                 dof_tmp = dof(ix);
%                 p = psi{dofNo}(x(keep,:));
%                 s_tmp = dof_tmp.*p;
%                 
%                 if nnz(keep) == size(double(s))
%                     s = s + s_tmp;
%                 else
%                     s(keep) = s(keep) + s_tmp;
%                 end                
                
                s(keep) = s(keep) + dof(ix).*psi{dofNo}(x(keep,:));
                
%                 ix = disc.getDofIx(state, dofNo, cells);
%                 s = s + dof(ix).*psi{dofNo}(x).*keep;
            end
            
        end
        
        %-----------------------------------------------------------------%
        function state = getCellSaturation(disc, state)
            
            [x, w, nq, ii, jj, cellNo] = makeCellIntegrator(disc.G, (1:disc.G.cells.num)', max(disc.degree), 'volume');
            W = sparse(ii, jj, w);

            x = disc.transformCoords(x, cellNo);
            
            sdof = state.sdof;
            nPh = size(sdof,2);
            s = zeros(disc.G.cells.num, nPh);
            for phNo = 1:nPh
                s(:,phNo) = (W*disc.evaluateSaturation(x, cellNo, sdof(:,phNo), state))./disc.G.cells.volumes;
            end
            
            state.s = s;
            
        end
        
        %-----------------------------------------------------------------%
        function I = cellInt(disc, integrand, cells, state)
        
            G    = disc.G;
            psi  = disc.basis.psi;
            grad_psi = disc.basis.grad_psi;
            nDof = state.nDof;
            nDofMax = disc.basis.nDof;
            
%             [x, w, nq, ii, jj, cellNo] = makeCellIntegrator(G, cells, max(disc.degree+1), 'volume');
%             W = sparse(ii, jj, w);
            
            
%             I = integrand(zeros(numel(cells)*nDofMax,disc.dim), ones(numel(cells)*nDofMax, 1), 1);
            
            I = double2ADI(zeros(sum(nDof),1), ...
                           integrand(zeros(sum(nDof),disc.dim), ones(sum(nDof), 1), 1, ones(1, disc.dim) ));
            
            for dofNo = 1:nDofMax
                
                keepCells = nDof(cells) >= dofNo;
                
                if any(keepCells)
                
                    ix = disc.getDofIx(state, dofNo, cells(keepCells));

                    [x, w, nq, ii, jj, cellNo] = makeCellIntegrator(G, cells(keepCells), max(disc.degree+1), 'volume');
                    W = sparse(ii, jj, w);
                    [x, ~, scaling] = disc.transformCoords(x, cellNo);

                    p  = psi{dofNo}(x);
                    gp = grad_psi{dofNo}(x).*scaling;
                    
%                     s  = disc.evaluateSaturation(x, cells(keepCells) , sdof , state );
%                     s0 = disc.evaluateSaturation(x, cells(keepCells)), sdof0, state0);
                    
                    I(ix) = W*integrand(x, cellNo, p, gp);
                    
                elseif numel(cells) == disc.G.cells.num
                    
                    warning('No cells with %d dofs', dofNo);
                    
                end
                
            end
            
            I = disc.trimValues(I);
            
        end
        
        %-----------------------------------------------------------------%
        function I = cellIntDiv(disc, integrand, cells, state)
            
            G        = disc.G;
            grad_psi = disc.basis.grad_psi;
            nDof     = state.nDof;
            nDofMax  = disc.basis.nDof;
            
%             [x, w, nq, ii, jj, cellNo] = makeCellIntegrator(G, cells, max(disc.degree+1));
%             W = sparse(ii, jj, w);
            
            I = double2ADI(zeros(sum(nDof),1), ...
                           integrand(zeros(sum(nDof),disc.dim), ones(sum(nDof), 1), 1));
            for dofNo = 1:nDofMax
                
                keepCells = nDof(cells) >= dofNo;
                
                if any(keepCells)
                    
                    ix = disc.getDofIx(state, dofNo, cells(keepCells));

                    [x, w, nq, ii, jj, cellNo] = makeCellIntegrator(G, cells(keepCells), max(disc.degree+1), 'volume');
                    W = sparse(ii, jj, w);
                    [x, ~, scaling] = disc.transformCoords(x, cellNo);

                    gp = grad_psi{dofNo}(x).*scaling;
                    I(ix) = W*integrand(x, cellNo, gp);
                    
                elseif numel(cells) == disc.G.cells.num
                    
                    warning('No cells with %d dofs', dofNo);
                    
                end
                    
            end
            
            I = disc.trimValues(I);
            
        end
        
        %-----------------------------------------------------------------%
        function I = faceIntDiv(disc, integrand, cells, upc, state)
            
            G       = disc.G;
            psi     = disc.basis.psi;
            nDof    = state.nDof;
            nDofMax = disc.basis.nDof;

%             [x, w, nq, ii, jj, cellNo, faceNo] = makeFaceIntegrator(G, cells, max(disc.degree+1));
%             W = sparse(ii, jj, w);

            upCells_v = G.faces.neighbors(:,2);
            intf      = find(disc.internalConn);
            upCells_v(intf(upc)) = disc.N(upc,1);
%             upCells_v = upCells_v(faceNo);    
%             upCells_G = upCells_v;
            
%             [x_c, ~, ~] = disc.transformCoords(x, cellNo);
%             [x_v, ~, ~] = disc.transformCoords(x, upCells_v);
%             [x_G, ~, ~] = disc.transformCoords(x, upCells_G);
            
            x0       = zeros(sum(nDof), disc.dim);
            [c0, f0] = deal(ones(sum(nDof), 1));
            I = double2ADI(zeros(sum(nDof),1), ...
                           integrand(x0, x0, x0, c0, c0, c0, f0, 1));
            for dofNo = 1:nDofMax
                
                keepCells = nDof(cells) >= dofNo;
                
                if any(keepCells)
                    
                    ix = disc.getDofIx(state, dofNo, cells(keepCells)');

%                     [x, w, nq, ii, jj, cellNo, faceNo] = makeFaceIntegrator(G, cells(keepCells), max(disc.degree+1));
                    [x, w, nq, ii, jj, cellNo, faceNo] = makeCellIntegrator(G, cells(keepCells), max(disc.degree+1), 'surface');
                    W = sparse(ii, jj, w);

                    upCells_vtmp = upCells_v(faceNo);
                    upCells_G = upCells_vtmp;

                    [x_c, ~, ~] = disc.transformCoords(x, cellNo);
                    [x_v, ~, ~] = disc.transformCoords(x, upCells_vtmp);
                    [x_G, ~, ~] = disc.transformCoords(x, upCells_G);

                    p = psi{dofNo}(x_c);
                    I(ix) = W*(integrand(x_c, x_v, x_G, cellNo, upCells_vtmp, upCells_G, faceNo, p));
                    
                elseif numel(cells) == disc.G.cells.num
                    
                    warning('No cells with %d dofs', dofNo);
                    
                end
                
            end
            
            I = disc.trimValues(I);
            
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
            
            jump = accumarray(cells, jump > 0.2) > 0;
            
            
        end
        
        %-----------------------------------------------------------------%
        function state = limiter(disc, state)
            
            G = disc.G;
            nDofMax = disc.basis.nDof;

            [smin, smax] = disc.getMinMaxSaturation(state);
            
            jump = disc.getInterfaceJumps(state);
            
            tol = 1e-4;
            over  = smax > 1 + tol;
            under = smin < 0 - tol;
            
            outside = over | under;
            bad     = outside | jump;
%             bad = jump;
            
            state.outside = outside;
            state.jump = jump;

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
                
                
%                 cells_o = find(over);
%                 cells_u = find(under);
%                 
%                 % Reduce to first-order
%                 ix = disc.getDofIx((1 + G.griddim + 1):nDof, cells);
%                 sdof(ix) = 0;
%                 
%                 % Reduce linear dofs so that saturaion is within [0,1]
% %                 ix0 = disc.getDofIx(1                , cells);
%                 for dNo = 1:G.griddim
%                     
%                     ix0 = disc.getDofIx(1, cells_u);
%                     ix  = disc.getDofIx(dNo+1, cells_u);
%                     sdof(ix) = sign(sdof(ix)).*min(abs(sdof(ix0)), abs(sdof(ix)));
%                     
%                     ix0 = disc.getDofIx(1, cells_o);
%                     ix  = disc.getDofIx(dNo+1, cells_o);
%                     sdof(ix) = sign(sdof(ix)).*(1 - sdof(ix0));
%                     
%                 end
%                 
%                 state.sdof(:,1) = sdof;
%                 state.sdof(:,2) = -sdof;
%                 
%                 ix0 = disc.getDofIx(1, (1:G.cells.num)');
%                 state.sdof(ix0,2) = 1 - state.sdof(ix0,1);
%                 
%                 [smin, smax] = disc.getMinMaxSaturation(state);
%                 
%                 s = disc.getCellSaturation(state);
                
            end
            
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
                
            s = disc.evaluateSaturation(xk, cellNo, state.sdof);
            
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
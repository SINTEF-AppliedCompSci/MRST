classdef DGDiscretization < WENODiscretization
    
    properties
        degree
        basis
        limiter
    end
    
    methods
        %-----------------------------------------------------------------%
        function disc = DGDiscretization(model, dim, varargin)
            
            disc = disc@WENODiscretization(model, dim);
            
            disc.degree  = 1;
            disc.basis   = 'legendre';
            disc.limiter = 'tvb';
            disc         = merge_options(disc, varargin{:});
            
%             disc.basis = DGBasisFunctions(disc.G, disc.degree);
            
            disc.basis   = dgBasis(disc.G, disc.degree, disc.basis);
            disc.limiter = dgLimiter(disc     , disc.limiter);
            
        end
        
        %-----------------------------------------------------------------%
        function [xhat, translation, scaling] = transformCoords(disc, x, cells)
            
            G           = disc.G;
            translation = -G.cells.centroids(cells,:);
            if isfield(G.cells, 'dx')
                scaling = 1./(G.cells.dx(cells,:)/2);
            else
                scaling     = 1./(G.cells.diameters(cells)/(2*sqrt(G.griddim)));
            end
            xhat        = (x + translation).*scaling;
               
        end
        
        %-----------------------------------------------------------------%
        function s = evaluateSaturation(disc, x, cells, dof)
            
%             if nargin < 5
%                 transCoord = true;
%             end
%             
%             if ~transCoord
%                 x = disc.transformCoords(x, cells);
%             end
            
            psi  = disc.basis.psi;
            nDof = disc.basis.nDof;
            
            s = 0;
            for dofNo = 1:nDof
                ix = (cells-1)*nDof + dofNo;
                s = s + dof(ix).*psi{dofNo}(x);
            end
            
        end
        
        %-----------------------------------------------------------------%
        function state = getCellSaturation(disc, state)
            
            ix = 1:disc.basis.nDof:disc.G.cells.num*disc.basis.nDof;
            
            cells = (1:disc.G.cells.num)';
            x = disc.G.cells.centroids;
            x = disc.transformCoords(x, cells);
            
            nPh = size(state.sdof,2);
            s = zeros(disc.G.cells.num, nPh);
            for phNo = 1:nPh
                s(:, phNo) = evaluateSaturation(disc, x, cells, state.sdof);
            end
            
            state.s = s;
            
        end
        
        %-----------------------------------------------------------------%
        function I = cellInt(disc, integrand, cells)
        
            G    = disc.G;
            psi  = disc.basis.psi;
            nDof = disc.basis.nDof;
            
            [x, w, nq, ii, jj, cellNo] = makeCellIntegrator(G, cells, disc.degree*G.griddim);
            W = sparse(ii, jj, w);
            
%             [x, cellNo, W] = cellBasisIntegrator(disc);
            
            [x, ~, ~] = disc.transformCoords(x, cellNo);
            
            I = integrand(repmat([0,0], numel(cells).*nDof, 1), ones(numel(cells)*nDof, 1), 1);
            for dofNo = 1:nDof
                ix = (1:nDof:numel(cells)*nDof) + dofNo - 1;
                p = psi{dofNo}(x);
%                 ix = (1:nDof:G.cells.num*nDof) + dofNo - 1;
                I(ix) = W*integrand(x, cellNo, p);
            end
            
        end
        
        %-----------------------------------------------------------------%
        function I = cellIntDiv(disc, integrand)
            
            G        = disc.G;
            grad_psi = disc.basis.grad_psi;
            nDof     = disc.basis.nDof;
            
            [x, cellNo, W] = cellBasisIntegrator(disc);
            
            [x, ~, scaling] = disc.transformCoords(x, cellNo);
            
            I = integrand(repmat([0,0], G.cells.num.*nDof, 1), ones(G.cells.num*nDof, 1), [1,1]);
            for dofNo = 1:nDof
                ix = (1:nDof:G.cells.num*nDof) + dofNo - 1;
                gp = grad_psi{dofNo}(x).*scaling;
                I(ix) = W*(integrand(x, cellNo, gp));
            end
            
        end
        
        %-----------------------------------------------------------------%
        function I = faceIntDiv(disc, integrand, upc)
            
            G        = disc.G;
            psi      = disc.basis.psi;
            nDof     = disc.basis.nDof;

            [x, cellNo, faceNo, W] = faceBasisIntegrator(disc);
            upCells_v = G.faces.neighbors(:,2);
            intf = find(disc.internalConn);
            upCells_v(intf(upc)) = disc.N(upc,1);
            upCells_v = upCells_v(faceNo);    
            upCells_G = upCells_v;
            
            [x_c, ~, ~] = disc.transformCoords(x, cellNo);
            [x_v, ~, ~] = disc.transformCoords(x, upCells_v);
            [x_G, ~, ~] = disc.transformCoords(x, upCells_G);
            
            x0 = repmat([0,0], G.cells.num.*nDof, 1);
            [c0, f0] = deal(ones(G.cells.num*nDof, 1));
            I = integrand(x0, x0, x0, c0, c0, c0, f0, 1);
            for dofNo = 1:nDof
                ix = (1:nDof:G.cells.num*nDof) + dofNo - 1;
                p = psi{dofNo}(x_c);
                I(ix) = W*(integrand(x_c, x_v, x_G, cellNo, upCells_v, upCells_G, faceNo, p));
            end
            
        end        
        
    end
        
end
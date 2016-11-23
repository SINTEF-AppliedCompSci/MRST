function [state, r] = conserveFlux(state, G, rock,  varargin)

    %   Merge input parameters
    opt = struct('bc'              , []       , ...
                 'src'             , []       , ...
                 'faceWeights'     , 'permWeighted', ...
                 'tol'             , 1e-14         );

    opt = merge_options(opt, varargin{:});
    
    %   Identify neumann boundary faces.
    neu = false(G.faces.num, 1);
    bf = boundaryFaces(G);
    neu(bf) = true;
    if ~isempty(opt.bc)
        neu(opt.bc.face(strcmp(opt.bc.type, 'pressure'))) = false;
    end
    
    %   Face orientations
    f = G.cells.faces(:,1);
    ncf = diff(G.cells.facePos);
    fSgn = 1 - 2*(G.faces.neighbors(f,1) ~= rldecode((1:G.cells.num)', ncf,1));
    
    %   Set rhs.
    rhs = zeros(G.cells.num,1);
    if ~isempty(opt.src)
        rhs(opt.src.cell) = opt.src.rate;
    end
    
    %   Matrix P sums faces per cell.
    [ii, jj] = blockDiagIndex(ones(G.cells.num,1), ncf);
    P = sparse(ii, jj, 1);
    
    %   Calculate residulas.
    r = rhs-P*(state.flux(f).*fSgn);

    %   If above tolerance, apply postprocessing.
    den = rhs;
    if nnz(rhs) == 0
        den = state.flux;
    end
    if norm(r)/norm(den) < opt.tol
        warning('Flux already conservative. No need for postprocessing.');
    else

        %   Calculate weights omega
        omega = computeFaceWeights(G, rock, opt);

        %   Build matrix systems as if all half-face normals have
        %   orientation out of the cell.
        
        %   Off-diagonal element (i,j) equals minus the face area of common
        %   face for cell i and j.
        c = G.faces.neighbors(f,:);
        c(fSgn == -1,:) = c(fSgn == -1, 2:-1:1);
        nz = c(:,2) ~= 0;
        B = sparse(c(nz,1), c(nz,2), -G.faces.areas(f(nz)).*omega(nz), ...
                                                 G.cells.num, G.cells.num);
        
        %   Diagonal element (i,i) equals sum of face areas for cell i.
        I = ones(numel(f),1);
        I(neu(f)) = 0;
        [ii, jj] = blockDiagIndex(ones(G.cells.num,1), ncf);
        ca = sparse(ii,jj,I)*(G.faces.areas(f).*omega);
        B = B + spdiags(ca, 0, G.cells.num, G.cells.num);
        
        if nnz(neu) == numel(bf)
            B(1,1) = 2*B(1,1);
        end

        beta = B\r;
        beta = rldecode(beta, ncf,1);
        beta(neu(f)) = 0;

        beta = sparse(f, 1:numel(f), I)*(omega.*beta.*fSgn);
        flux = state.flux + G.faces.areas.*beta;

        r = rhs-P*(flux(f).*fSgn);
        if norm(r)/norm(den) > opt.tol;
            warning('Could not construct conservative flux field');
        end

        state.flux = flux;
        
    end

end

function omega = computeFaceWeights(G, rock, opt)

    if strcmp(opt.faceWeights, 'permWeighted')
        
        f = G.cells.faces(:,1);
        ncf = diff(G.cells.facePos);
        fSgn = 1 - 2*(G.faces.neighbors(f,1) ~= rldecode((1:G.cells.num)', ncf,1));
        
        K = permTensor(rock, G.griddim)';
        [ii, jj] = blockDiagIndex(G.griddim*ones(G.cells.num,1));
        K = sparse(ii, jj, K(:));
        
        fn = bsxfun(@times,G.faces.normals(f,:), fSgn./G.faces.areas(f))';
        [ii, jj] = blockDiagIndex(G.griddim*ones(G.cells.num,1), ncf);
        fnMat = sparse(ii,jj, fn(:));
        delta = fnMat'*K*fnMat;
        delta = spdiags(delta, 0);

        ii = f;
        jj = (1:numel(f))';
        omega = 1./(sparse(ii, jj,1)*delta);

        for i = 1:G.faces.num
            d = delta(f == i);
            omega(i) = omega(i)*(numel(d)*prod(d));
        end
        
        omega = omega(f);
        
    elseif strcmp(opt.faceWeights, 'tpf')
        
        omega = computeTrans(G, rock)./G.faces.areas(G.cells.faces(:,1));
        
    elseif strcmp(opt.faceWeights, 'ones')
        omega = ones(numel(G.cells.faces(:,1)),1);
    end
    
end
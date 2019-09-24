function [points,weights,nPts] = fitMoments3(x, basis, moments, varargin)

    opt = struct('reduce', true, 'tol', 1e-10, 'chunkSize', 5);
    opt = merge_options(opt, varargin{:});

    psi  = basis.psi;
    nDof = basis.nDof;
    n0 = size(x,1);
    
    ne = numel(moments{1});

    nPts   = zeros(ne,1);
    [points, weights] = deal([]);

    x0     = x;
    reduce = opt.reduce;
    nc = floor(ne/opt.chunkSize);
    ePos   = round(linspace(0, ne, nc))+1;
    for i = 1:numel(ePos)-1
        
        elements = ePos(i):ePos(i+1)-1;
        ne_loc = numel(elements);
        
        fprintf('Compressing quadrature for element %d to %d of %d ... ', ...
                                         min(elements), max(elements), ne);
        tic;
        
        m = cellfun(@(m)m(elements), moments, 'UniformOutput', false);
        tol = 100*eps(mean(m{1}));
        m = vertcat(m{:});
        M = nan(numel(m),1);
        for j = 1:nDof
            M((1:nDof:nDof*ne_loc)+j-1) = m((ne_loc*(j-1) + (1:ne_loc)));
        end
        
        reduced = true;
        x = x0;
        n = n0;
        k = n0;
        w = zeros(n,1);
        
        while k > 0 && reduced 
        
            % Matrix of basis functions evalauted at current quadrature points            
            P = computeBasisMatrix(psi, x, ne_loc, n, nDof);
            
            s = sum(full(P(1:nDof, 1:n)).^2,1);
            if reduce
                [~, ix] = sort(s);
            else
                ix = 1;
            end
            reduced = false;
            xPrev = x;
            wPrev = w;
            nPrev = n;

            for pNo = 1:numel(ix)
                x(ix(pNo),:) = [];
                n = size(x,1);
                
                P = computeBasisMatrix(psi, x, ne_loc, n, nDof);
                
                if 0
                    w = P'*((P*P')\M);
                    flag = all(w > 0);
                else
                    [w, ~, residual, flag] = lsqnonneg(P, M);
                    flag = flag && all(abs(residual)<tol);
                end
                if flag > 0
                    k       = k-1;
                    reduced = true;
                    break
                else
                    x = xPrev;
                    w = wPrev;
                    n = nPrev;
                end    
            end
            
            if ~opt.reduce
                break
            end
            
        end
        
        nPts(elements) = n;
        points  = [points;repmat(x,ne_loc,1)];
        weights = [weights;w];
        
        time = toc;
        fprintf('Compressed from %d to %d points in %f second\n', n0, n, time);
              
    end

    
end

function P = computeBasisMatrix(psi, x, ne, n, nDof)
    p = cell2mat(cellfun(@(p) p(x), psi, 'unif', false));
    [ii, jj] = blockDiagIndex(n*ones(ne,1), nDof*ones(ne,1));
    P = sparse(jj, ii, repmat(p, ne,1));
end
function [x,w,nPts] = fitMoments(x, basis, moments, num)

    psi  = basis.psi;
    nDof = basis.nDof;
    nPts = size(x,1);
    
    k = nPts;
    reduced = true;
    w = [];
    
    opts = optimoptions('lsqlin');
    opts.ConstraintTolerance = 1e-12;
    opts.OptimalityTolerance = 1e-6;
    opts.MaxIterations       = 100;
    
    while k > 0 && reduced
    
        % Matrix of basis functions evalauted at current quadrature points
        p = cell2mat(cellfun(@(p) p(x), psi, 'unif', false));
        s = sum(reshape(p.^2, nPts, nDof), 2);
        [~, ix] = sort(s);
        reduced = false;
        xPrev = x;
        wPrev = w;
        
        for pNo = 1:numel(ix)
            x(ix(pNo),:) = [];
            nPts = size(x,1);
            p = cell2mat(cellfun(@(p) p(x), psi, 'unif', false));
            [ii, jj] = blockDiagIndex(nPts*ones(num,1), nDof*ones(num,1));
            P = sparse(jj, ii, repmat(p, num,1));

            I = speye(nPts*num); o = zeros(nPts*num,1); r = Inf(nPts*num,1);
            [w, ~, ~, flag] = lsqlin(I, o, I, r, P  , moments, [], [], [], opts);
            
            if flag > 0
                k       = k-1;
                reduced = true;
                break
            else
                x = xPrev;
                w = wPrev;
            end    
        end
        
    end
    
    nPts = size(x,1);
    
end
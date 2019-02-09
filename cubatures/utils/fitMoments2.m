function [points,weights,nPts] = fitMoments2(x, basis, moments, varargin)

    opt = struct('equal', false, 'reduce', true, 'tol', 1e-12);
    opt = merge_options(opt, varargin{:});

    psi  = basis.psi;
    nDof = basis.nDof;
    n = size(x,1);
    
%     opts = optimoptions('lsqlin');
%     opts.ConstraintTolerance = 1e-12;
%     opts.OptimalityTolerance = 1e-6;
%     opts.MaxIterations       = 100;
%     opts.Display             = 'off';
    opts = optimset('lsqnonneg');
    opts = optimset('TolFun', 1e-12);
    
    if ~opt.equal
        nElements = numel(moments{1});
    else
        nElements = 1;
    end

    nPts   = zeros(nElements,1);
    [points, weights] = deal(cell(nElements,1));

    x0     = x;
    n0     = n;
    reduce = opt.reduce;
    parfor cNo = 1:nElements
        
        fprintf('Compressing quadrature for element %d of %d ... ', cNo, nElements);
        tic;
        
        M = cellfun(@(m)m(cNo), moments);
        vol = M(1);
        M = M./vol;
        reduced = true;
        x = x0;
        n = n0;
        k = n0;
        w = zeros(n,1);
        
        while k > 0 && reduced 
        
         % Matrix of basis functions evalauted at current quadrature points
            P = reshape(cell2mat(cellfun(@(p) p(x), psi, 'unif', false)), n, nDof)';
            s = sum(P.^2,1);
            if reduce
                [~, ix] = sort(s);
            else
                ix = [];
            end
            reduced = false;
            xPrev = x;
            wPrev = w;
            nPrev = n;

            for pNo = 1:numel(ix)
                x(ix(pNo),:) = [];
                n = size(x,1);
                P = reshape(cell2mat(cellfun(@(p) p(x), psi, 'unif', false)), n, nDof)';

                I = eye(n); o = zeros(n,1); r = Inf(n,1);
%                 [w, ~, ~, flag] = lsqlin(I, o, I, r, P  , M, [], [], [], opts);
%                 [w, ~, ~, flag] = lsqlin(I, o, I, r, P  , M, o, r, [], opts);
                [w, ~, residual, flag] = lsqnonneg(P, M);
%                 tol = 1e-15;
                flag = flag && all(abs(residual)<opt.tol);
%                 PP = P*P';
%                 w = P'*(PP\M);
%                 tol = 0.01;
%                 flag = all(isfinite(w)) && rank(P) == nDof && abs(rcond(PP)) > tol;
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
        
        nPts(cNo)    = n;
        points{cNo}  = x;
        weights{cNo} = w.*vol;
        
        time = toc;
        fprintf('Compressed from %d to %d points in %f second\n', n0, n, time);
              
    end

    
end
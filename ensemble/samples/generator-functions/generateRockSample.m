function sample = generateRockSample(N, varargin)
    opt = struct('avg_perm'         , 100*milli*darcy, ...
                 'std_perm'         , 0.5            , ...
                 'avg_poro'         , 0.5            , ...
                 'std_poro'         , 0.1            , ...
                 'min_poro'         , 1e-3           , ...
                 'seed'             , []             , ...
                 'covarianceFun'    , []             , ...
                 'correlationLength', 0.3            , ...
                 'smooth'           , false          , ...
                 'toVector'         , false          );
    opt = merge_options(opt, varargin{:});

    fun = setCovarianceFun(opt);
    
    if ~isempty(opt.seed)
        rng(opt.seed);
    end
    
    if numel(N) - nnz(N == 1) == 1
        p = GaussianProcess1D(N, fun);
    else
        p = GaussianProcessND(N, fun);
    end
    
    p = p - mean(p(:));
    p = p./std(p(:));
    
    perm = 10.^(log10(opt.avg_perm) + p.*opt.std_perm);
    poro = opt.avg_poro + p.*opt.std_poro;
    poro = max(poro, opt.min_poro);
    
    if opt.toVector
        perm = perm(:); poro = poro(:);
    end
    
    sample = struct('perm', perm, 'poro', poro);
end

function fun = setCovarianceFun(opt)
    fun = opt.covarianceFun;
    if isempty(opt.covarianceFun)
        if ~opt.smooth
            fun = @(x) exp(-sqrt(sum(x.^2,2))/opt.correlationLength);
        else
            fun = @(x) exp(-sum(x.^2,2)/opt.correlationLength);
        end
    end
end
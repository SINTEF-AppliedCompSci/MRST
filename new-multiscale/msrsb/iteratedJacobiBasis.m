function [I, interactionMap, res] = iteratedJacobiBasis(A, CG, varargin)
    opt = struct('enforceCenter', false, ...
                 'iterations',    round(10*CG.parent.cells.num/CG.cells.num), ...
                 'incrementTol',  1e-3, ...
                 'autostop',      true, ...
                 'verbose',       mrstVerbose(), ...
                 'useConstant',   true, ...
                 'omega',         2/3, ...
                 'interpolator',  [], ...
                 'limitSupport',  true ...
                 );
    opt = merge_options(opt, varargin{:});
    centers = mapCenters(CG);
    CG.cells.centers = centers;
    
    
    A = A - diag(sum(A, 2));
    
    if ~isfield(CG.cells, 'interaction')
        CG = storeInteractionRegion(CG);
    end
    
    lens = cellfun(@numel, CG.cells.interaction);
    blocks = rldecode((1:CG.cells.num)', lens);
    interaction = vertcat(CG.cells.interaction{:});
    interactionMap = sparse(interaction, blocks, ones(size(interaction)), CG.parent.cells.num, CG.cells.num);

    G = CG.parent;
    
    isNode = false(G.cells.num, 1);
    if isempty(opt.interpolator)
        if opt.useConstant
            I_0 = controlVolumeRestriction(CG.partition)';
        else
            I_0 = linearInterpolator(CG);
        end
    else
        I_0 = opt.interpolator;
    end
    
    if opt.enforceCenter
        isNode(centers) = true;

        rhs = -A(~isNode, isNode)*speye(numel(centers));
        A = A(~isNode, ~isNode);
        interactionMap = interactionMap(~isNode, :);
    else
        rhs = 0*I_0;
    end
    I = I_0(~isNode, :);
    
    D = diag(A);
    n = numel(D);
    D_inv = spdiags(1./D, 0, n, n);
    
    res = nan(min(opt.iterations, 1000), 1);
    w = opt.omega;
    i = 1;
    while i < opt.iterations
        def = D_inv*(A*I - rhs);
        if opt.limitSupport
             def = removeOutOfBounds(def, interactionMap);
        end
        nval = full(sum(abs(def), 1));
        if size(res, 1) < i
            res = [res; nan*res];
        end
        
        res(i) = norm(nval, 2);
        if i > 1 && (opt.autostop && res(i) > res(i-1) || abs(res(i)-res(i-1)) < opt.incrementTol)
            dispif(opt.verbose, 'Finished after %d iterations\n', i);
            break
        end
        dispif(opt.verbose, '%d of %d iterations in basis\n', i, opt.iterations)
        
        I = I - w*def;
        I = bsxfun(@rdivide, I, sum(I, 2));
        i = i + 1;
    end
    I_0(~isNode, :) = I;
    I = I_0;
end

function increment = removeOutOfBounds(increment, interactionMap)

    increment = increment.*interactionMap;
end
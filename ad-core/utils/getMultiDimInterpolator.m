function fn = getMultiDimInterpolator(x, Y)
% Get a multidimensional interpolator (with support for ADI varibles)
    sz = size(Y);
    assert(iscell(x));
    if numel(x) > 1
        assert(numel(sz) == numel(x));
        assert(all(cellfun(@numel, x) == sz));
    end
    
    nvar = numel(x);
    dY = cell(nvar, 1);
    
    T = griddedInterpolant(x, Y, 'linear', 'linear');
    for i = 1:nvar
        dyi = diff(Y, 1, i);
        dxi = diff(x{i});
        

        
        % delta y / delta x
        if nvar > 1
            % Permute to allow for use of bsxfun
            p = 1:nvar;
            p(i) = 1;
            p(1) = i;
            dyidx = bsxfun(@rdivide, permute(dyi, p), dxi);
            % Permute back
            dyidx = ipermute(dyidx, p);
        else
            dyidx = dyi./dxi;
        end

        % Evaluate using midpoints
        xi = x;
        xi{i} = xi{i}(1:end-1) + dxi/2;
        dY{i} = griddedInterpolant(xi, dyidx, 'nearest', 'nearest');
    end
    
    fn = @(varargin) interpTableND(T, dY, varargin{:});
end

function Ye = interpTableND(Y, dY, varargin)
    Yq = varargin;
    nvar = numel(dY);
    assert(numel(Yq) == nvar);
    
    isad = cellfun(@(x) isa(x, 'ADI'), Yq);
    
    if any(isad)
        % Get a sample variable
        ad = Yq(isad);
        ad = ad{1};
        
        % Create double as input for the function evaluation
        Yqd = cellfun(@double, Yq, 'UniformOutput', false);
        Ye = Y(Yqd{:});
        % Cast to ADI
        Ye = double2ADI(Ye, ad);
        for i = 1:nvar
            if isad(i)
                % If it is AD, compute derivative using chain rule
                dydx = dY{i}(Yqd{:});
                ix = (1:numel(dydx))';
                d = sparse(ix, ix, dydx);
                for j = 1:numel(Ye.jac)
                    Ye.jac{j} = Ye.jac{j} + d*Yq{i}.jac{j};
                end
            end
        end
    else
        % Just use interpolator directly
        Ye = Y(Yq{:});
    end
end
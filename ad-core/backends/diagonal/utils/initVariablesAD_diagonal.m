function varargout = initVariablesAD_diagonal(varargin)
% Diagonal AD initializer
    assert(nargout == nargin || nargout == nargin - 1);
    
    values = varargin(1:nargout);
    types = ones(nargin, 1);
    numVals = cellfun(@numel, values);
    if nargin > nargout
        types = varargin{end};
    else
        for i = 2:nargin
            if numVals(i) ~= numVals(i-1)
                types(i) = types(i-1) + 1;
            else
                types(i) = types(i-1);
            end
        end
    end
    offset = zeros(max(types)+1, 1);
    ctr = 2;
    offset(1) = 1;
    for i = 2:numel(types)
        if types(i) ~= types(i-1)
            offset(ctr, 1) = i;
            ctr = ctr+1;
        end
    end
    offset(end) = numel(types) + 1;
    
    ntypes = max(types);

    n = nargout;
    varargout = cell([1, n]);
    for i = 1:n
        zerojac = cell(1, ntypes);
        for j = 1:ntypes
            type = j;
            sub = types == type;
            nv = numVals(sub);
            nval = nv(1);
            dim = [nval, nnz(sub)];
            d = zeros(numVals(i), 0);
            zerojac{j} = DiagonalJacobian(d, dim, zeros(numVals(i), 1));
        end
        varargout{i} = NewAD(varargin{i}, zerojac);
        varargout{i}.numVars = numVals';
        varargout{i}.offsets = offset;
    end
    for type = 1:ntypes
        sub = find(types == type);
        nv = numVals(sub);
        nval = nv(1);
        nsub = numel(sub);
        
        assert(all(nv == nval));
        for i = 1:nsub
            s = sub(i);
            tmp = zeros(nval, nsub);
            tmp(:, i) = 1;
            
            dim = [nval, nsub];
            varargout{s}.jac{type} = DiagonalJacobian(tmp, dim);
        end
    end
end
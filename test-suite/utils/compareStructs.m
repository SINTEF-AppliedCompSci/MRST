function out = compareStructs(in1, in2, varargin)
    
    opt = struct('fun'     , @abs , ...
                 'relative', false, ...
                 'omit'    , {{}} , ...
                 'includeStructs', true);
    opt = merge_options(opt, varargin{:});
    
    if numel(in1) > 1
        assert(numel(in1) == numel(in2))
        out = [];
        for i = 1:numel(in1)
            o   = compareStructs(in1(i), in2(i), varargin{:});
            out = [out; o]; %#ok
            return;
        end
    end
    
    assert(isstruct(in1) && isstruct(in2), 'Both inputs must be structs');
    
    keep1 = getValidFields(in1, opt);
    keep2 = getValidFields(in2, opt);
    
    names1 = fieldnames(in1); names1 = names1(keep1);
    names2 = fieldnames(in2); names2 = names2(keep2);
    
    names   = intersect(names1, names2);
    names12 = union(names1, names2);
    
    unhandled = setdiff(names12, names);
    if ~isempty(unhandled)
        warning('Found %d fields not common to both structs', numel(unhandled));
    end
    
    out = cell2struct(num2cell(nan(numel(names12),1)), names12);
    
    for i = 1:numel(names)
        name = names{i};
        v1 = in1.(name);
        v2 = in2.(name);
        if isstruct(v1) && opt.includeStructs()
            assert(isstruct(v2))
            out.(name) = compareStructs(v1, v2, varargin{:});
        else
            if iscell(v1)
                assert(iscell(v2))
                v = cellfun(@(v1,v2) compareValues(v1, v2, opt), ...
                                           v1, v2, 'UniformOutput', false);
            else
                v = compareValues(v1, v2, opt);
            end
            out.(name) = v;
        end
    end
    
end

%-------------------------------------------------------------------------%
function keep = getValidFields(in, opt)
    v = struct2cell(in);
    n = numel(v);
    keep = false(n, 1);
    keepfun = @(v) isnumeric(v) || islogical(v);
    if opt.includeStructs
        keepfun = @(v) keepfun(v) || isstruct(v);
    end
    for i = 1:n
        if iscell(v{i})
            keep(i) = all(cellfun(@(v) keepfun(v), v{i}));
        else
            keep(i) = keepfun(v{i});
        end 
    end
    if ~isempty(opt.omit)
        keep = keep & cellfun(@(n) ~any(strcmpi(n, opt.omit)), fieldnames(in));
    end
end

%-------------------------------------------------------------------------%
function v = compareValues(v1, v2, opt)
    v = opt.fun(v1 - v2);
    if opt.relative
        v = v./max(max(opt.fun(v1), opt.fun(v2)), 1e-10);
    end
end

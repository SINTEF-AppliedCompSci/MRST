function [G, names, category, keep] = buildStateFunctionDigraph(C, names, category, varargin)
    opt = struct('Start', {{}}, ...
                 'Stop', {{}}, ...
                 'Center', {{}}, ...
                 'FilterNames', {names});
    opt = merge_options(opt, varargin{:});

    n = numel(names);
    G = digraph(C, names);
    keep = true(n, 1);
    
    hasStart = ~isempty(opt.Start);
    hasStop = ~isempty(opt.Stop);
    hasCenter = ~isempty(opt.Center);
    
    if hasStop || hasCenter
        D = distances(G);
    end
    
    if hasStart || hasCenter
        Gt = digraph(C', names);
        Dr = distances(Gt);
    end
    
    if hasStart
        keep = keep & filter(Dr, opt.FilterNames, opt.Start);
    end
    
    if hasStop
        keep = keep & filter(D, opt.FilterNames, opt.Stop);
    end
    
    if hasCenter
        keep = keep & (filter(D,  opt.FilterNames, opt.Center) | ...
                       filter(Dr, opt.FilterNames, opt.Center));
    end
    
    G = G.subgraph(keep);
    names = names(keep);
    category = category(keep);
end

function keep = filter(D, names, targets)
    keep = false(numel(names), 1);
    if ischar(targets)
        targets = {targets};
    end
    for i = 1:numel(targets)
        name = targets{i};
        pos = strcmpi(names, name);
        keep = keep | isfinite(D(:, pos));
    end
end
function [ws, T] = convertJutulWellSols(wells)
    T = wells.time;
    names = setdiff(fieldnames(wells), 'time');
    n = numel(T);
    nw = numel(names);
    ws = cell(n, 1);
    
    first = names{1};
    has_water = isfield(wells.(first), 'wrat');
    has_oil = isfield(wells.(first), 'orat');
    has_gas = isfield(wells.(first), 'grat');

    default = struct('name', 'dummy', 'bhp', 0, 'status', true);
    if has_water; default.qWs = 0; end
    if has_oil; default.qOs = 0; end
    if has_gas; default.qGs = 0; end

    ws0 = repmat(default, 1, nw);
    for i = 1:n
        w = ws0;
        for j = 1:nw
            name = names{j};
            W = wells.(name);
            w(j).name = name;
            w(j).bhp = W.bhp(i);
            if has_water
               w(j).qWs = W.wrat(i);
            end
            if has_oil
               w(j).qOs = W.orat(i);
            end
            if has_gas
               w(j).qGs = W.grat(i);
            end
        end
        ws{i} = w;
    end
end

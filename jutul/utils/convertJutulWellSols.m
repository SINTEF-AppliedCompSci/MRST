function [ws, T] = convertJutulWellSols(wells)
    T = wells.TIME;
    names = setdiff(fieldnames(wells), 'TIME');
    n = numel(T);
    nw = numel(names);
    ws = cell(n, 1);
    
    first = names{1};
    has_water = isfield(wells.(first), 'WRAT');
    has_oil = isfield(wells.(first), 'ORAT');
    has_gas = isfield(wells.(first), 'GRAT');

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
            w(j).bhp = W.BHP(i);
            if has_water
               w(j).qWs = W.WRAT(i);
            end
            if has_oil
               w(j).qOs = W.ORAT(i);
            end
            if has_gas
               w(j).qGs = W.GRAT(i);
            end
        end
        ws{i} = w;
    end
end

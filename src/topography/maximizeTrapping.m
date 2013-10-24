function trees = maximizeTrapping(G, varargin)
    opt = struct('res', [], 'n', inf);
    opt = merge_options(opt, varargin{:});
    
    if isempty(opt.res)
        opt.res = trapAnalysis(G);
    end
    
    v = volumes(G, opt.res, 1:max(opt.res.traps));

    % All root nodes - traps that are downstream from every other trap
    rn = rootNotes(opt.res.trap_adj');
    Nt = min(opt.n,numel(rn));
    
    trees = struct('root', [], 'traps', [], 'value', cellfun(@(x) 0, cell(1,Nt), 'UniformOutput', false));
    % Initialize
    i = 1;
    
    downstream = arrayfun( @(ind) find(dfs(opt.res.trap_adj, ind) > -1), rn, ...
                                                  'UniformOutput', false);
    rn_vols = cellfun(@(x) sum(v(x)), downstream);
    while i <= Nt
        [vol, index] = max(rn_vols);
        
        ds = downstream{index};
        trees(i).value = vol;
        trees(i).traps = ds;
        trees(i).root = rn(index);
        v(ds) = 0;
        
        rn_vols = cellfun(@(x) sum(v(x)), downstream);
        i = i + 1;
    end
end

function v = volumes(G, res, trap)
% Find volume of a subset of traps
    v = zeros(1,numel(trap));
    for i = 1:numel(trap)
        ind = res.traps == trap(i);
        z = G.cells.z(ind);
        fill_z = res.trap_z(trap(i));
        v(i) = sum(max(eps, G.cells.volumes(ind).*(fill_z - z)));
    end
end

function c = rootNotes(A)
    c = find(sum(A,2) == 0);
end

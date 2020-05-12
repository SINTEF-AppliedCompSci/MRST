function T = TransNTPFA(G, u, OSflux, varargin)
% Set up using homogeneous Neumann bcs
    
    opt = struct('avgmpfa', false);
    opt = merge_options(opt, varargin{:});

    dispif(mrstVerbose, 'TransNTPFA, AvgMPFA=%d\n', opt.avgmpfa);
    
    T = cell(2, 1);
    
    m = G.faces.num;
    n = G.cells.num;
    
    internal = 1:m;
    internal(~all(G.faces.neighbors ~= 0, 2)) = [];

    tii = cell(2, 1);
    tsp = cell(2, 1);
    r = cell(2, 1);
    mu = cell(2, 1);
    
    for j = 1:2
        tend = zeros(m, 1);
        tii{j} = zeros(m, 2);
        nij = zeros(m, 1);
        
        % Sweep for sparsity set up
        for i = internal
            t = OSflux{i, j};
            nij(i) = max(0, size(t, 1) - 3);
            tii{j}(i, 1:2) = [t(1, 2), t(2, 2)];
        end
        
        % Set up sparsity pattern and values
        ii = zeros(sum(nij), 1);
        jj = zeros(sum(nij), 1);
        vv = zeros(sum(nij), 1);
        s = [1; cumsum(nij) + 1];
        
        % Don't loop over zero rows
        nzrows = (1:m);
        nzrows(nij == 0) = [];
        
        for i = nzrows
            t = OSflux{i, j};
            idx = s(i):(s(i+1) - 1);
            ii(idx) = i;
            jj(idx) = t(3:end-1, 1);
            vv(idx) = t(3:end-1, 2);
            tend(i) = t(end, 2);
        end
        
        tsp{j} = sparse(ii, jj, vv, m, n);
        
        r{j} = tsp{j} * u + tend;
        mu{j} = 0 * r{j} + 0.5;
        T{j} = 0 * r{j};
    end
    
    if ~opt.avgmpfa
        epstol = 1e-12 * max(full(max(tsp{1}, [], 2)), full(max(tsp{2}, [], 2)));
        for j = 1:2
            ir = abs(r{j}) <= epstol;
            r{j}(ir) = 0;
        end
        jj = abs(r{1} + r{2}) > epstol;
        mu{1}(jj) = r{2}(jj) ./ (r{1}(jj) + r{2}(jj)); 
        mu{2}(jj) = ones(sum(jj), 1) - mu{1}(jj);
        assert(all(mu{1} >= 0.0))
        assert(all(mu{2} >= 0.0))
        assert(all(mu{1} <= 1.0))
        assert(all(mu{2} <= 1.0))
    end

    T{1}(internal) = mu{1}(internal) .* tii{1}(internal, 1) + mu{2}(internal) .* tii{2}(internal, 2);
    T{2}(internal) = mu{1}(internal) .* tii{1}(internal, 2) + mu{2}(internal) .* tii{2}(internal, 1);
end


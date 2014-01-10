function dg = cull_dual_partition(g, dg, varargin)
   opt = struct('verbose',      mrstVerbose,   ...
                'cg',    [], ...
                'CullEdges', true,...
                'SystemMatrix', [],...
                'Iter', nan, ...
                'Cull', true, ...
                'FixEdges', true ...
                 );

   opt = merge_options(opt, varargin{:});

    % Find indices of the edge cells - cells not in the "inner" partition
    cells = 1:g.cells.num;


    tic()
    if ~opt.CullEdges
        if isempty(opt.cg)
            error('Coarse grid must be supplied to avoid culling edges')
        end
        cg = opt.cg;
        edge_blocks = cg.faces.neighbors(any(cg.faces.neighbors==0,2),:);
        edge_blocks = edge_blocks(edge_blocks ~=0);
        % Create a mask which ignores fine cells if their parent coarse
        % block is on the edge of the domain
        inside_mask = ~ismember(cg.partition, edge_blocks);
    else
        inside_mask = true(g.cells.num,1);
    end

    tmp = true(g.cells.num,1);
    tmp(dg.ii) = false;
    tmp(dg.lineedge) = false;
    e_cells = cells(tmp & inside_mask);
    clear tmp;
    if isempty(opt.SystemMatrix)
        % Generate graph for the system
        inner = g.faces.neighbors(:,1) ~= 0 & g.faces.neighbors(:,2) ~= 0;
        tmp = double(g.faces.neighbors(inner,:));
        A1 = sparse(vertcat(tmp(:,1), tmp(:,2), (1:g.cells.num)'),...
                    vertcat(tmp(:,2), tmp(:,1), (1:g.cells.num)'),1);
    else
        A1 = opt.SystemMatrix ~= 0;
    end


    % We want to select the neighbors of two or more neighbors to
    % handle cartesian grids
    A = A1 | (A1*A1==2);
    fprintf('Constructed system graph in %2.2fs\n', toc())
    % Create a copy of the graph, but remove entries corresponding to nodes
    % which are not in the "inner" partition
    A_inner = A1;

    tic();
    A_inner(e_cells,:) = 0;
    fprintf('Removed edges in %2.2fs\n', toc())
    tic();
    comp = components(A_inner);
    fprintf('Found components in %2.2fs\n', toc())
    % Set category of edge cells to zero. These are not strongly connected
    % components in anything but the trivial sense, and will thusly have
    % been assigned unique indices we don't need
    comp(e_cells) = 0;
    % Set inner nodes to category -1 to ensure that they are never
    % connected to the inner partition, which would violate mass balance
    % when doing msfv calculations on the grid
%     comp(dg.nn) = -1;

    tic();
    ii_new = zeros(1, numel(e_cells));
    [compi, compj] = find(A);
    % Never recategorize nodes directly connected to inner nodes
    [~, skiplist] = find(A1(dg.nn,:));
    skiplist = [skiplist; dg.lineedge(:)];
    iter = 0;

    e_cells = e_cells(randperm(numel(e_cells)))
%     e_cells = e_cells(end:-1:1);
    while opt.Cull
        changed = zeros(numel(e_cells),1);
        for c = 1:numel(e_cells)
            c_i = e_cells(c);
            if any(skiplist == c_i)
                continue
            end
            % Find the categories of all connected neighbors by using the
            % original graph to find the neighbor structure
            comp_neighbors = comp(compj(compi==c_i));
            pool = unique(comp_neighbors(comp_neighbors ~= 0));
            if numel(pool) == 1 % && ~any(c_i == dg.nn)
                ii_new(c) = c_i;
                comp(c_i) = pool;
                changed(c) = 1;
%                 clf; plotCellData(g, comp); plotGrid(g, 'FaceColor', 'None'); colormap colorcube; pause();
            end
        end
        e_cells = e_cells(changed==0);
        numel(e_cells)
        iter = iter + 1;
        fprintf('Iteration %d culled %d cells\n', iter, sum(changed));
        if all(changed==0) || opt.Iter == iter
            break;
        end
    end

    if opt.FixEdges
        tmp = false(g.cells.num,1);
        tmp(dg.nn) = true;
        tmp(dg.lineedge) = true;
        A_outer = A1;
        A_outer(~tmp,:) = 0;
        comp_outer = components(A_outer);
        bad_lineedge = find(~ismember(comp_outer(dg.lineedge), comp_outer(dg.nn)));
        fprintf('Removed %d bad linear edges...\n', numel(bad_lineedge));
        dg.lineedge = setdiff(dg.lineedge, bad_lineedge);
    end


    ii_new = ii_new(ii_new ~= 0);
    dg.ii = [dg.ii ii_new];

    dg.ee = dg.ee(~ismember(dg.ee, ii_new) );
    dg.lineedge = dg.lineedge(~ismember(dg.lineedge, ii_new));

    fprintf('Culled %d cells in %2.2fs (Converged after %d cycles) \n', numel(ii_new), toc(), iter)
end

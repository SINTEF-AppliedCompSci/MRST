function dg = dualPartitionGraph(cg, varargin)

   opt = struct(...
       'targets', 1:cg.cells.num,...
       'virtual', true,...
       'cache', false...
   );
   % cache off may not give correct linear edges
   opt = merge_options(opt, varargin{:});

    dg.ee = []; dg.ii = zeros(cg.cells.num,1); dg.nn = []; dg.lineedge = [];
    neighbors = cell(cg.cells.num,1);
    for i = 1:cg.cells.num
        [inner ~] = find(cg.faces.neighbors == i);
        neigh = unique(cg.faces.neighbors(inner,:));
        neighbors{i} = neigh(neigh>0)';
    end
    if opt.virtual
        [cg, neighbors] = addVirtualCG(cg, neighbors, 10);
    else
        cg.is_virtual = zeros(cg.cells.num,1);
    end

    for i = 1:cg.cells.num
        tmp = neighbors{i};
        % Because addvirtual and this function have a bit different
        % approach to wether a node is it's own neighbor
        neighbors{i} = tmp(tmp~=i);
    end

    NC = cg.cells.num;

    % system matrix
    g = cg.parent;
    inner = g.faces.neighbors(:,1) ~= 0 & g.faces.neighbors(:,2) ~= 0;
    tmp = double(g.faces.neighbors(inner,:));
    A = sparse(vertcat(tmp(:,1), tmp(:,2), (1:g.cells.num)'),...
                vertcat(tmp(:,2), tmp(:,1), (1:g.cells.num)'),1);


    centers = zeros(cg.cells.num,1);
    for i = 1:cg.cells.num
        pt = cg.cells.centroids(i,:);
        dist = bsxfun(@minus, g.cells.centroids, pt);
        dist = sqrt(sum(dist.^2,2));
        [val ind] = min(dist);
        % Find the index of the nearest fine cell to the coarse block
        % centroid
        centers(i) = ind;
    end

    tmp = [];
    for i = 1:cg.cells.num
        fprintf('%d of %d\n', i, cg.cells.num)
        a = centers(i);
        [v pred] = shortest_paths(A,a);
        neigh = neighbors{i};
        for j = 1:numel(neigh)
            b = centers(neigh(j));
            if cg.is_virtual(i) && cg.is_virtual(neigh(j))
                continue
            end
            p = path_from_pred(pred, b);
            tmp = [tmp p];
%             clf
%             plotGrid(g, p)
%             plotGrid(g, neigh(j), centers(i))
        end
    end
    c = count_coarse_edges(cg, A);
    dg.lineedge = tmp';

    t = false(g.cells.num, 1);

    bnd = g.faces.neighbors(sum(g.faces.neighbors == 0,2)==1,:) ;
    bnd = bnd(bnd>0);
    bnd = bnd(:);

    t(bnd(bnd>0)) = true;


    counts = accumarray(bnd,1);

%     dg.ii = find(c==3);
%
    dg.ii = find(c == 4 |...
                (c == 3 & t) |...
                (c == 2 & counts > 1)| ...
                counts > 2);

    dg.nn = centers(~cg.is_virtual);
    dg.nn = unique(dg.nn);

    dg.ee = [];
    dg = cullDual(g, dg, A);

    dg.lineedge = setdiff(dg.lineedge, dg.ii);
    dg.lineedge = setdiff(dg.lineedge, dg.nn);
    dg.ii = setdiff(dg.ii, dg.nn);
end

function c = count_coarse_edges(CG, A)
    % c(i) is the number of coarse blocks fine cell i is once removed from
    Cn = CG.cells.num;
    cn = CG.parent.cells.num;
    cnts = false(cn, Cn);
    for i = 1:Cn
        tmp = false(cn, 1);
        tmp(CG.partition == i) = true;
        % Find neighboring connections
        cnts(:,i) = A*tmp > 0;
    end
    c = sum(cnts,2);
end

function i_cent = snap_to_cell(CG, coarseindex)
    i_cent  = CG.cells.centroids(coarseindex,:);
    if ~CG.is_virtual(coarseindex)
        [gi li] = findCenter(CG, i_cent, find(CG.partition==coarseindex));
        % Use not the centroid, but the centroid of the fine cell nearest
        % to the centroid
        i_cent = CG.parent.cells.centroids(gi,:);
    end
end

function n = common_neighbors(ci, cj, neighbors)
    %ci scalar, cj scalar
    n = intersect(neighbors{ci}, neighbors{cj});
end

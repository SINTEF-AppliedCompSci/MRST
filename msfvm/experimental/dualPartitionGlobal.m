function DG = dualPartitionGlobal(CG, varargin)

   opt = struct(...
       'targets', 1:CG.cells.num,...
       'virtual', true,...
       'virtualmod', 1, ...
       'cache', false...
   );
   % cache off may not give correct linear edges
   opt = merge_options(opt, varargin{:});

    DG.ee = []; DG.ii = zeros(CG.cells.num,1); DG.nn = []; DG.lineedge = [];
    neighbors = cell(CG.cells.num,1);
    for i = 1:CG.cells.num
        [inner ~] = find(CG.faces.neighbors == i);
        neigh = unique(CG.faces.neighbors(inner,:));
        neighbors{i} = neigh(neigh>0)';
    end
    if opt.virtual
        [CG, neighbors] = addVirtualCG(CG, neighbors, opt.virtualmod);
    else
        CG.is_virtual = zeros(CG.cells.num,1);
    end

    for i = 1:CG.cells.num
        tmp = neighbors{i};
        % Because addvirtual and this function have a bit different
        % approach to wether a node is it's own neighbor
        neighbors{i} = tmp(tmp~=i);
    end

    NC = CG.cells.num;
    processed = zeros(10*NC,3);
    next_p = 1;

    nt = numel(opt.targets);
    for ii = 1:nt
        i = opt.targets(ii);
        fprintf('%d of %d \n',i, nt);
        neigh = neighbors{i};
        i_cent = snap_to_cell(CG, i);
        n_cents = zeros(numel(neigh),3);
        %% Find centroids and align them to nearest cell centroid to get better grid orientation
        for ni = 1:numel(neigh)
           n_cents(ni,:) = snap_to_cell(CG, neigh(ni));
        end
        ee = cell(numel(neigh),1);
        le = cell(numel(neigh),1);
        for ni = 1:numel(neigh)
            %% Find common neighbors
            cn = common_neighbors(i, neigh(ni), neighbors);
            if isempty(cn) || 1
                % Look at neighbors of neighbors common to the first
                for nj = 1:numel(neigh)
                    if neigh(ni) == neigh(nj)
                       continue
                    else
                        % If the two neighbors of the current coarse blocks
                        % have a common neighbor, it should be considered
                        % (For example north and south block in Cart grids)
                        overlap = intersect(neighbors{neigh(ni)}, neighbors{neigh(nj)});
                        overlap = overlap(overlap ~= i);
                        overlap = overlap(~ismember(overlap, neighbors{i}));
                        if ~isempty(overlap)
                            cn = [cn neigh(nj)]; %#ok
                        end
                    end
                end
                cn = cn(cn ~= i);
            end
            c1 = n_cents(ni,:);
            N_cn = numel(cn);
            ee_new = cell(N_cn,1);
            le_new = cell(N_cn,1);
            for nj = 1:N_cn
                c_tuple = sort([i, neigh(ni), cn(nj)]);
                %% Caching stuff - DRY
                if opt.cache &&...
                   any(processed(:,1) == c_tuple(1) &...
                       processed(:,2) == c_tuple(2) &...
                       processed(:,3) == c_tuple(3))
                    continue
                else
                    processed(next_p, :) = c_tuple;
                    next_p = next_p + 1;
                    if next_p>size(processed,1)
                        processed = vertcat(processed, zeros(10*NC, 3)); %#ok, intentional expansion
                    end
                end
                if all(CG.is_virtual(c_tuple))
                   continue
                end
                %% Actually categorizing the nodes
                % For all pairs, categorize the nodes using the category
                % function based on the centroids
                c2 = snap_to_cell(CG, cn(nj));
                f = create_category_function(i_cent, c1, c2);

%                 blockInd = find(CG.partition == i | CG.partition == neigh(ni) | CG.partition == cn(nj));
                %TODO: flytt ut "finn noder"- biten for nabolag
                blockInd = find(CG.partition == i | CG.partition == neigh(ni) | ismember(CG.partition, cn));

                [ee_new{nj} le_new{nj}] = categorize(blockInd, CG, f);
%                 [ee_new{nj} le_tmp] = categorize(blockInd, CG, f);
%                 if any(neighbors{neigh(ni)}==neigh(nj))
%                    le_new{nj} = le_tmp;
%                 end
            end

            if ~opt.cache
                for nj = 1:N_cn
                    for nk = nj+1:N_cn
                        isect = intersect(le_new{nj}, le_new{nk});
                        if ~isempty(isect)
                            for nl = nk+1:N_cn
                                isect = intersect(isect, le_new{nl});
                                if ~isempty(isect)
                                    for nm = nl+1:N_cn
                                        DG.lineedge = [DG.lineedge intersect(isect, le_new{nm})];
                                    end
                                end
                            end
                        end
                    end
                end
            end
            ee{ni} = horzcat(ee_new{:});
            le{ni} = horzcat(le_new{:});
        end

        DG.ee =[DG.ee horzcat(ee{:})];
%         DG.lineedge =[DG.lineedge horzcat(le{:})];
        if ~CG.is_virtual(i)
            DG.nn(i) = findCenter(CG, i_cent, find(ismember(CG.partition,i)));
        end
    end
    DG.nn = DG.nn(DG.nn ~= 0);
    DG.ee = unique(DG.ee);
    DG.ii = setdiff(1:CG.parent.cells.num, [DG.ee DG.nn]);
    DG.ee = setdiff(DG.ee, DG.nn);
    if opt.cache
        DG = dual_add_linear_edges(CG.parent, DG);
    else
        DG.lineedge = unique(setdiff(DG.lineedge, DG.nn));
    end
end


function Plane = createPlane(center, node1, node2)
    N = cross(node1-center, node2-center);
    Plane = @(pt) N*(pt - repmat(center,size(pt,1),1))';
end

function [ee_new le_new] = categorize(blockInd, cg, f)
    N = numel(blockInd);
    [p n] = cellnode(cg.parent, blockInd);
    coords = cg.parent.nodes.coords(n,:);
    orient = f(coords);

    le = cell(N,1);
    ee = cell(N,1);

    real_orient = real(orient);
    imag_orient = imag(orient);
    for c = 1:N
        nodePos = p(c):(p(c+1)-1);
        oI =     imag_orient(nodePos);
        % If there is an imaginary component, we are inside the correct
        % part of the domain for the orientation to be valid.
        if any(oI)
            oR = real_orient(nodePos);
            if ~(all(oR>0)|| all(oR<0))
               ee{c} = blockInd(c);
               if ~all(oI);
                   le{c} = blockInd(c);
               end
            end

        end
    end
    ee_new = horzcat(ee{:});
    le_new = horzcat(le{:});
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


function categoryf = create_category_function(center, p1, p2)
    % Normal vector to the plane
    N = cross(p1-center, p2-center);
    N = N./norm(N);

    plane   = createPlane(center, p1, p2);
    pts = {center, p1, p2};
    d = 3;
    p = cell(d,1);
    for ind = 1:d
        ii = mod(ind+1,3)+1;
        jj = mod(ind+3,3)+1;
        % Create a plane which passes through two of the current points,
        % and is orthogonal to the current plane to divide the domain
        p{ind} = orthoplane(pts{ind}, pts{jj}, pts{ii}, N);
    end
    categoryf = @(pt) plane(pt) + (p{1}(pt) & p{2}(pt) & p{3}(pt))*1i;
end

function f = orthoplane(p_current, p_other, p_ortho, N)
    tmp = createPlane(p_current, p_other, p_current + N);
    sign_current = sign(tmp(p_ortho));
    f = @(pt) sign(tmp(pt)) == sign_current;
end

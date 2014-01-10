function dual = dualPartition(cg, varargin)
% Create coarse dual partition/grid
%
% SYNOPSIS:
%    Creates a dual grid from a coarse grid
% DESCRIPTION:
%    Uses a variety of partitioning algorithms to create a dual grid for
%    different coarse grids
%
% REQUIRED PARAMETERS:
%    cg -- coarse grid as defined by ex. partitionUI
%
% OPTIONAL PARAMETERS:
%    Scheme -- Partitioning scheme to be used. Different schemes include:
%              - 'plane' uses a plane to divide the coarse grid. Only works
%              for right angled coarse grid blocks
%              - 'polyplane' uses two planes and different selection algorithm.
%              Generally better at complex grids.
%              - 'polyplane2' same as polyplane, but targets cell centers
%              instead of points in space. A bit better for most grids.
%              Generally this is your best bet for unstructured grids
%              - 'orthoplane' planar using centroids and directions along
%              the principal axes. Not recommended for complex grids.
%              - For logically cartesian grids, consider using
%              partitionUIdual which is very fast.
%    StoreEdge -- Store edges
%
% RETURNS:
%    dual -- Categorization of the nodes
% NOTE:
%    There is no surefire way to get a working dual grid for a given grid.
%    Experiment with different coarse grids and algorithms to get a good
%    result.
    opt = struct('Verbose', mrstVerbose, 'Scheme', 'polyplane', 'StoreEdge', false);
    opt = merge_options(opt, varargin{:});

    N = cg.cells.num;
    dual.nn = zeros(1,N);
    dual.ii = [];
    dual.ee = [];
    dual.edges = cell(N,cg.griddim);
    dual.lineedge = [];

    if opt.Verbose; h = waitbar(0,sprintf('0 of %d coarse volumes processed', N)); end;
    for i = 1:N
        %find corresponding nodes
        blockInd = find(cg.partition == i);
        facenum = numel(cg.cells.facePos(i): cg.cells.facePos(i+1)-1);
        if facenum < 4
            % Sometimes impossible coarse blocks exist in unstructured
            % grids. We ignore them.
            continue
        end
        if strcmp(opt.Scheme, 'plane') || facenum < 6
            func = partitionPlanes(dual, cg, i, blockInd);
            fprintf('%d faces and %d cells in coarse block\n', facenum, numel(blockInd));
        elseif strcmp(opt.Scheme, 'orthoplane')
            c = cg.cells.centroids(i,:);
            func = cell(3,1);
            func{1} = createPlane(c, c + [1 0 0], c + [0 0 1]);
            func{2} = createPlane(c, c + [1 0 0], c + [0 1 0]);
            func{3} = createPlane(c, c + [0 1 0], c + [0 0 1]);
        else
            func = partitionPoly(cg, i, blockInd, opt.Scheme);
        end
        dual.nn(i) = findCenter(cg, cg.cells.centroids(i,:), blockInd);

        % Find face positions for all cells
        c_start = cg.parent.cells.facePos(blockInd);
        c_stop  = cg.parent.cells.facePos(blockInd+1)-1;
        % Number of faces for each cell
        c_sizes = c_stop-c_start + 1;
        fa = cg.parent.cells.faces(mcolon(c_start, c_stop));

        f_start = cg.parent.faces.nodePos(fa);
        f_stop = cg.parent.faces.nodePos(fa+1)-1;
        % Number of nodes for each face
        f_sizes = f_stop-f_start + 1;
        nodeIndices = cg.parent.faces.nodes(mcolon(f_start , f_stop));

        % index into the faces
        f_index = 1;
        % index into nodes (which is what we want)
        n_index = 1;
        nodeInd = cell(numel(blockInd),1);
        for c = 1:numel(blockInd)
            numberofnodes = sum(f_sizes(f_index:(f_index+c_sizes(c)-1)));
            nodeInd{c} = n_index:(n_index+numberofnodes-1);
            n_index = n_index + numberofnodes;
            f_index = f_index + c_sizes(c);

        end
        %%
        coords =  cg.parent.nodes.coords(nodeIndices,:);
        orient = zeros(size(coords,1),cg.griddim);
        for pp = 1:cg.griddim
            % Save the orientation of all points
           orient(:,pp) = func{pp}(coords);
        end

        for c = 1:numel(blockInd)
            nodePos = nodeInd{c};
            for pp = 1:cg.griddim
                orientation = orient(nodePos, pp);
                if ~(all(orientation>0)|| all(orientation<0))
%                 if abs(sum(sign(orientation))) < numel(orientation)
%                 if sum(orientation > 0) ~= numel(orientation) && sum(orientation < 0) ~= numel(orientation)
                    dual.edges{i,pp} =  [dual.edges{i,pp} blockInd(c)];
                    if ~opt.StoreEdge
                        % We will not need the intersections and can break
                        % prematurely
                        break;
                    end
                end
            end
        end
        if cg.griddim == 3 && opt.StoreEdge
            dual.lineedge = [ dual.lineedge intersect(dual.edges{i,1}, dual.edges{i,2}), ...
                                 intersect(dual.edges{i,2}, dual.edges{i,3}), ...
                                 intersect(dual.edges{i,1}, dual.edges{i,3})];
        end
        if opt.Verbose;waitbar(i/N,h,sprintf('%d of %d coarse volumes processed',i, N)); end;
    end
    dual.ee = horzcat(dual.edges{:});
    dual.nn = dual.nn(dual.nn~=0);

    if opt.Verbose; close(h); end;
    dual.ee = unique(dual.ee);
    %remove the center nodes and zero element from neigbhour array
    dual.ee = setdiff(dual.ee, [0 dual.nn]);
    if cg.griddim == 3 && opt.StoreEdge
        if ~ismember(0, dual.lineedge)
            dual.lineedge = unique(dual.lineedge);
        else
            dual.lineedge = unique(dual.lineedge(dual.lineedge > 0));
        end
    end
    dual.ii = setdiff(1:cg.parent.cells.num, [dual.ee dual.nn]);
end


function nb = findNeighbourFaces(coarseindex, CG)
    %finds three neighbouring faces which are all each others neighbours to
    %be used for plane generation

    %get the indices of the faces for this coarse cell
    faces = CG.cells.faces( CG.cells.facePos(coarseindex): ...
                            (CG.cells.facePos(coarseindex+1) - 1),1);

    %the idea is to find a single pivot point and from this find the faces
    %which touch this pivot point such that these three faces must be
    %neighbours. For a somewhat non-degenerate system, this can be used to
    %find intersecting planes for partitioning the dual grid around the
    %centroid

    Nf = numel(faces);
    %select the first face arbitrarily
    nb = [faces(1)];
    nodelist = getFineNodeCoords(faces(1), CG);

    loop = 0;
    while(1)
        %for ind = 2:Nf
        for ind = setdiff(1:Nf, nb(1))
            facenodes = getFineNodeCoords(faces(ind), CG);
            isect = intersect(facenodes, nodelist);
            if isect
                nodelist = isect;
                nb = unique([nb faces(ind)]);
                if numel(nb) == 3
                    break;
                end
            end
        end
        loop = loop + 1;
        if numel(nb) == CG.griddim
            break;
        end
        if loop > 3
            %if for some reason the loop doesn't converge, restart with
            %another starting point...
            t = floor(Nf*rand() + 1);
            nb = [faces(t)];
            nodelist = getFineNodeCoords(faces(t), CG);
        end
    end
end

function N = getFineNodeCoords(faceindex, CG)
    %Given a coarse face index, find all the corresponding fine node
    %positions
    finefaces = CG.faces.fconn(CG.faces.connPos(faceindex) : CG.faces.connPos(faceindex + 1) - 1);
    %get the indices of the positions
    nodepos = CG.parent.faces.nodePos(sort([finefaces' finefaces'+1]));
    %subtract one from every second element to get correct endpoints
    nodeindices = [];
    for i = 1:numel(nodepos)/2
        nodeindices = [nodeindices nodepos(i*2 - 1):(nodepos(i*2)-1)];
    end
    %find the node positions
    N = unique(CG.parent.faces.nodes(nodeindices));
end

function Plane = createPlane(center, node1, node2)
    %create a normal to the plane using the cross product definition
    N = cross(node1-center, node2-center);
    %a plane is defined by N dot (plane point - any point) = 0. We can
    %exploit this to get a function defining which side of a plane some
    %random point is
    Plane = @(pt) N*(pt - repmat(center,size(pt,1),1))';
end

function Plane = createPolyPlane(center, top, bottom, start, stop, type)
    % Find a vector which is inbetween the top and start vectors
    Nc = (top - center + start - center)/2;
    % Use this vector along with the center to define a plane which divides
    % the domain
    div = @(pt) Nc*(pt - repmat(center,numel(pt)/3,1))';
    % Check the orientation of the plane and points and define divider
    % accordingly
    if div((top + start + center)/3) > 0
        Divider = @(pt) div(pt) > 0;
    else
        Divider = @(pt) div(pt) < 0;
    end
    switch regexprep(type,'\d','')
        case 'polynomial'
            % Create Vandermonde matrix for linear polynomials for both the planes
            tmp = [center; top; bottom; start; stop];
            x = tmp(:,1);
            y = tmp(:,2);
            z = tmp(:,3);
            V = [x y z ...
                 x.*x y.*y z.*z ...
                 x.*y x.*z y.*z ones(5,1)];
            % Penrose inverse for least squares fitting
            c1 = V\ones(5,1);

            Plane = @(pt)   ([pt(:,1), pt(:,2), pt(:,3), ...
                              pt(:,1).*pt(:,1), pt(:,2).*pt(:,2), pt(:,3).*pt(:,3), ...
                              pt(:,1).*pt(:,2), pt(:,1).*pt(:,3), pt(:,2).*pt(:,3), ...
                              ones(numel(pt)/3,1)]*c1)' - 1;

        case 'polyplane'
            P1 = createPlane(center, start, top);
            P2 = createPlane(center, stop, bottom);
            Plane = @(pt)  Divider(pt).*P1(pt) + ...
                          ~Divider(pt).*P2(pt);
        otherwise
            error('Error!')
    end
end

function planes = generatePlanes(center, node1, node2, node3)
    %three nodes and a center defines three planes
    if numel(center) == 3
        planes = cell(3,1);
        planes{1} = createPlane(center, node1, node2);
        planes{2} = createPlane(center, node1, node3);
        planes{3} = createPlane(center, node2, node3);
    else
        planes = cell(2,1);
        planes{1} = linefunc(node1, center);
        planes{2} = linefunc(node2, center);
    end
    return
end

function planes = generatePolyPlanes(center, a1, a2, b1, b2, top, bottom, type)
    %three nodes and a center defines three planes
    planes = cell(3,1);
    planes{1} = createPolyPlane(center, top, bottom, a1, a2, type);
    planes{2} = createPolyPlane(center, top, bottom, b1, b2, type);
    planes{3} = createPolyPlane(center, b1, b2, a1, a2, type);
    return
end

function l = linefunc(node1, node2)
    l = @(pt)(node2(:,1)-node1(:,1))*(pt(:,2)-node1(:,2))-(node2(:,2)-node1(:,2))*(pt(:,1)-node1(:,1));
    return
end

function planes = partitionPoly(cg, i, blockInd, type)
    f = cg.cells.faces(cg.cells.facePos(i):cg.cells.facePos(i+1)-1);
    [nsub, sub] = subFaces(cg.parent, cg);
    sub_ix = cumsum([0; nsub]);

    if strcmp(type, 'polyplane2')
        centAlt = 1;
    else
        centAlt = 0;
    end
    cent = cg.cells.centroids(i,:);
    if centAlt
        c = findCenter(cg, cent, blockInd);
        cent = cg.parent.cells.centroids(c,:);
    end
    cents = cg.faces.centroids(f',:);
    areas = cg.faces.areas(f,:);
    norms = cg.faces.normals(f',:);
    neighbors = cg.faces.neighbors(f,:);
    % Transform the normals so they are relative to the current block
    norms = norms.*repmat(2*(neighbors(:,1)==i) - 1, 1, 3);
    % Define top and bottom face by taking the largest and smallest z
    % vector
    [cM topFace] = max(norms(:,3));
    [cm underFace] = min(norms(:,3));
    topC = cents(topFace,:);
    underC = cents(underFace,:);
    % ensure that top and bottom face won't get picked
    areas([topFace underFace]) = 0;
    % Sort the faces by area and find the local indices of the four largest
    % faces to be used in pairs
    % Big faces are selected to increase the odds of hitting "real" faces
    % and not just degenerate areas. This will also hopefully make the
    % connections consistent across primal coarse blocks
    [sortedArea, sortIndex] = sort(areas,'descend');
    largestFaces = sortIndex(1:4);
    %find first face and pair it up with the face with the biggest diff in
    %angle
    a1 = largestFaces(1);
    dp = zeros(3,1);
    for j = 2:4
        dp(j-1) = acos(dot(norms(a1,1:2), norms(largestFaces(j), 1:2)));
    end
    [tmp1 mind] = max(dp);
    a2 = largestFaces(mind + 1);
    rest = setdiff(largestFaces, [a1 a2]);
    b1 = rest(1);
    b2 = rest(2);
    if centAlt
        t = @(face) faceToBlockCent(cg, cents(face,:), sub(sub_ix(f(face)) + 1 : sub_ix(f(face) + 1)));
        planes = generatePolyPlanes(cent, t(a1), t(a2), t(b1), t(b2), t(topFace), t(underFace), type);
    else
        planes = generatePolyPlanes(cent, cents(a1,:), cents(a2,:), cents(b1,:), cents(b2,:), topC, underC, type);
    end
end

function planes = partitionPlanes(dual, cg, i, blockInd)
    %% Partition using planes and three pivot points - fallback for edge coarse blocks
    %find the center nodes
    center = cg.cells.centroids(i,:);
    dual.nn(i) = findCenter(cg, center, blockInd);
    neighbours = findNeighbourFaces(i, cg);
    centroids = cg.faces.centroids(neighbours,:);
    if cg.griddim == 3
        planes = generatePlanes(center, centroids(1,:), centroids(2,:), centroids(3,:));
    else
        planes = generatePlanes(center, centroids(1,:), centroids(2,:), 0);
    end
    return
end

function newcent = faceToBlockCent(cg, facecent, faces)
    facecentroids = cg.parent.faces.centroids(faces,:);
    if cg.griddim == 3
        %find distance from coarse centroid to fine centroids
        centroidDist =  (facecentroids(:,1)-facecent(1)).^2 + ...
                        (facecentroids(:,2)-facecent(2)).^2 + ...
                        (facecentroids(:,3)-facecent(3)).^2;
    else
        centroidDist =  (facecentroids(:,1)-facecent(1)).^2 + ...
                        (facecentroids(:,2)-facecent(2)).^2;
    end
    [val least] = min(centroidDist);
    globalFaceIndex = faces(least);
    neighbors = setdiff(cg.parent.faces.neighbors(globalFaceIndex,:),0);
    cellsN = cg.parent.cells.centroids(neighbors,:);

    if numel(neighbors) == 2
        newcent = sum(cellsN)/2;
    else
        newcent = cellsN;
    end
end

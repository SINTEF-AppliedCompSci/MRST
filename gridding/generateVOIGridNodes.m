function [pSurfs, t, bdyID] = generateVOIGridNodes(GC, packed, WR, layerRf, opt)
% Generate 3D points of all surfaces of volume of interest (VOI) and 
% connectivity list corresponding to the 2D planar points
%
% SYNOPSIS:
%   [pSurfs, t, bdyID] = generateVOIGridNodes(GC, packed, WR, layerRf, opt)
%
% PARAMETERS:
%  GC           - CPG grid structure
%  packed       - Structure of VOI information, see 'allInfoOfVolume'
%  WR           - Structure of 2D WR points, see 'prepareWellRegionNodes2D'
%  layerRf      - Number of refinement layers in each VOI layer
%  opt          - Parameters for generating the unstructured grid, see
%                 'ReConstructToUnstructuredGrid'
%
% RETURNS:
%  pSurfs   - Points of all VOI grid surfaces
%  t        - Connectivity list
%  bdyID    - Indices of outer boundary nodes in points pSurfs
%
% SEE ALSO:
%  `VolumeOfInterest` `generateHWGridNodes`

    % Reference surface
    refSurf = ceil( length(packed.nodes)/2 );
    switch opt.gridType
        case 'triangular'
            generator = @triangularPts;
        case 'Voronoi'
            generator = @VoronoiPts;
        otherwise
            error([mfilename, ': Unknown grid type'])
    end
    
    % Get 2D points and conectivity list
    [p, t, bdyID] = generator(GC, packed, WR, refSurf, opt);
    
    % Get surface points
    pSurfs = getSurfacePoints(GC, packed, layerRf, refSurf, p, bdyID);
end

function [p, t, bdyID] = triangularPts(G, packed, WR, refSurf, opt)
    % Assign VOI data
    bnV = packed.bdyNodes{refSurf};
    
    % Assign WR data
    pW  = [vertcat(WR.points.cart); vertcat(WR.points.rad)];
    tW  = WR.connlist;
    bnW = WR.bdnodes;
    
    % Get boundary points
    pib = pW(bnW,:); % Inner boundary points
    pob = G.nodes.coords(bnV, [1,2]); % Outer boundary points
    
    % Generate basic points from 'DistMesh'
    [pdis, fd] = passToDistmesh(pib, pob, opt.multiplier, opt.maxIter);
    
    % Delaunay triangulation
    delTri  = delaunayTriangulation(pdis);
    
    % Remove cells indise the WR
    tol = 0.1;
    tV    = delTri.ConnectivityList;
    pVmid = ( pdis(tV(:,1),:) + pdis(tV(:,2),:) + pdis(tV(:,3),:) ) / 3; 
    in = fd(pVmid) < -tol;
    tV    = tV(in,:);
    
    % Add the WR points and connectivity list
    npib = size(pib, 1);
    npob = size(pob, 1);
    nIDIB = tV > npib;
    npWR  = size(pW,1);
    tV( nIDIB) = tV(nIDIB) - npib + npWR;
    tV(~nIDIB) = bnW( tV(~nIDIB) );
    tV = mat2cell(tV, ones(size(tV,1), 1));
    t  = [tW; tV];
    p  = [pW; pdis(npib+1 : end, :)];
    bdyID  = npWR + (1:npob)';
end

% -------------------------------------------------------------------------
function [p, t, bdyID] = VoronoiPts(G, packed, WR, refSurf, opt)
    % Get outer boundary points and auxiliary points of VOI
    % see 'gridding examples\connWRCartGrids'
    fV     = packed.faces{refSurf};
    bfV    = packed.bdyFaces{refSurf};
    boxfV  = packed.boxFaces{refSurf};
    bnV    = packed.bdyNodes{refSurf};
    pob    = G.nodes.coords(bnV, [1,2]);
    pob2   = G.faces.centroids(bfV, [1,2]);
    boxfV  = boxfV(~ismember(boxfV, fV));
    pauxV  = G.faces.centroids(boxfV, [1,2]); % box face centroids
    
    % Get inner boundary points and auxiliary points of WR
    % see 'gridding examples\connWRCartGrids'
    pW   = [vertcat(WR.points.cart); vertcat(WR.points.rad)];
    tW   = WR.connlist;
    bnW  = WR.bdnodes;
    pib  = pW(bnW,:);
    [pIn, pOut, R] = computeAuxPts(pW, bnW, 0.23);
    pib2  = pOut;
    pauxW = pIn;
    
    % Generate basic points from 'DistMesh'
    pdis = passToDistmesh(pib2, pob2, opt.multiplier, opt.maxIter, ...
        'pIBRadius', R);
    
    % Get Voronoi points and connectivity list
    pall = [pdis; pauxV; pauxW];
    [pVor, tVor] = voronoin(pall, {'Qbb','Qz'});
    
    % Clip the diagram
    fdI = @(p)dpoly(p, [pib; pib(1,:)]);
    fdO = @(p)dpoly(p, [pob; pob(1,:)]);
    fd  = @(p)ddiff(fdO(p), fdI(p));
    tol1 = 0.1;
    tol2 = 0.1;
    try
        [pVor, tVor] = clipDiagram(pVor, tVor, fd, tol1, tol2);
    catch
        [pVor, tVor] = clipDiagram2(pVor, tVor, fd, tol1, tol2);
    end
    
    % Add WR points, and map the connectivity list again
    D2 = euclideanDistance(pVor, pib);
    [IBID, ~] = find( bsxfun(@eq, D2, min(D2)) );
    map1 = find( ~ismember((1:size(pVor,1))', IBID) );
    pVor = pVor(map1, :);
    map1 = [(1:length(map1))', map1];
    map1(:,1) = map1(:,1) + size(pW,1);
    map2 = [bnW, IBID];
    tVor = cellfunUniOut(@(x)[map1(ismember(map1(:,2), x), 1)', ...
        map2(ismember(map2(:,2), x), 1)'], tVor);
    p = [pW; pVor];
    t = [tW; tVor];
    D3 = euclideanDistance(p, pob); 
    [bdyID, ~] = find( bsxfun(@eq, D3, min(D3)) );
    
    % Add empty cells
    [p, t] = addEmpCells(p, t, bnW);
end

% -------------------------------------------------------------------------
function players = getSurfacePoints(GC, packed, layerRf, refSurf, p, bdyID)
    % Generate the 3D points of all surfaces of VOI grid
    bnV  = packed.bdyNodes;
    nV   = packed.nodes;
    np   = size(p, 1);
    inID = find(~ismember((1:np)', bdyID));
    players = cell(length(nV), 1);
    
    % Get points of each VOI surfaces
    for k = 1 : length(nV)
        xi = p(inID, 1);
        yi = p(inID, 2);
        % Shift the inner points
        xcRef = mean( GC.nodes.coords(bnV{refSurf}, 1) );
        ycRef = mean( GC.nodes.coords(bnV{refSurf}, 2) );
        xcNow = mean( GC.nodes.coords(bnV{k}, 1) );
        ycNow = mean( GC.nodes.coords(bnV{k}, 2) );
        xi = xi + (xcNow - xcRef);
        yi = yi + (ycNow - ycRef);
        % Interpolate the Z-coords of inner points
        n = unique( cell2mat(nV{k}) );
        x = GC.nodes.coords(n, 1);
        y = GC.nodes.coords(n, 2);
        z = GC.nodes.coords(n, 3);
        zi = griddata(x, y, z, xi, yi);
        pfull = zeros(np,3);
        pfull(inID, [1, 2])    = [xi, yi];
        pfull(inID,  3    )    =  zi;
        pfull(bdyID, [1, 2, 3]) = GC.nodes.coords(bnV{k}, :);
        players{k} = pfull;
    end

    % Add points of refined surfaces
    for k = 1 : length(players)-1
        pUpp = players{k};
        pBot = players{k+1};
        scal = linspace(0, 1, layerRf(k)+1);
        scal = scal(2:end-1);
        prefine = arrayfunUniOut(@(xyz)(pBot(:,xyz)-pUpp(:,xyz))/(1-0)*(scal-0) ...
            + pUpp(:,xyz), (1:3)');
        prefine = cell2mat(prefine);
        prefine = arrayfunUniOut(@(col)reshape(prefine(:, col), [], 3), ...
            (1:size(prefine,2))');
        prefine = cell2mat(prefine);
        players{k} = [players{k}; prefine];
    end 
    players = cell2mat(players);
    nlayer  = size(players, 1) / np;
    players = arrayfunUniOut(@(L)players( (1:np) + (L-1)*np, : ), (1:nlayer)');
end

% -------------------------------------------------------------------------
function [pVor, tVor] = clipDiagram(pVor, tVor, fd, tol1, tol2)
    % Remove points outside the region
    in  = find(fd(pVor) < tol1 & all(~isinf(pVor), 2));
    map = [(1:length(in))', in];
    pVor = pVor(map(:,2), :);
    
    % Remove conflict points (point too close to each other)
    D = euclideanDistance(pVor, pVor);
    D = triu(D);
    
    [removed, reserved] = find(D < tol2);
    ii = removed < reserved; 
    removed = removed(ii);
    reserved = reserved(ii);
    map(removed,1) = map(reserved,1);
    idx  = find(~ismember((1:size(pVor,1))', removed));
    pVor = pVor(idx, :);
    map(:,1) = arrayfun(@(x)find(x == idx), map(:,1));
    
    % Map the connectivity list
    tVor = cellfunUniOut(@(x)unique( map(ismember(map(:,2), x), 1)' ), tVor);
    tVor = tVor( cellfun(@length, tVor) > 3 );
end

% -------------------------------------------------------------------------
function [pVor, tVor] = clipDiagram2(pVor, tVor, fd, tol1, tol2)
    % Remove points outside the region
    cCenter = cellfunUniOut(@(t)mean(pVor(t, :)), tVor);
    cCenter = cell2mat(cCenter);
    in = fd(cCenter) < tol1;
    t = tVor(in);
    n = cell2mat(t')';
    n = unique(n);
    p = pVor(n, :);
    t = cellfunUniOut(@(t)find(ismember(n, t)), t);
    t = sortPtsCounterClockWise(p, t);
    g = tessellationGrid(p, t);
    g = removeShortEdges(g, tol2);
    pVor = g.nodes.coords;
    tVor = arrayfunUniOut(@(c)gridCellNodes(g, c), (1:g.cells.num)');
end

% -------------------------------------------------------------------------
function [pIn, pOut, R] = computeAuxPts(p, bn, m0)
    pib = p(bn, :);
    pib = [pib; pib(1,:)];
    n   = size(bn,1);
    e2n = [(1:n)', [(2:n)';1]];
    % Compute the radius
    edges = [bn, [bn(2:end); bn(1)]];
    L  = p(edges(:,1),:) - p(edges(:,2),:);
    L  = sqrt(sum(L.^2,2));
    do = true;
    m  = m0;
    while do
        if m > 0.5
            throwError(L)
        end
        try
            m = m + 0.02;
            R = ( [L(1);L(1:n-1)] + [L(n);L(2:n)] ) * m;
            [pIn, pOut] = deal( zeros(size(edges,1), 2) );
            for i = 1 : size(edges,1)
                [x1, y1, x2, y2] = deal(p(edges(i,1),1), p(edges(i,1),2), ...
                    p(edges(i,2),1), p(edges(i,2),2));
                [r1, r2] = deal(R(e2n(i,1)), R(e2n(i,2)));
                pCross = circleCross(x1, y1, r1, x2, y2, r2);
                in = inpolygon(pCross(:,1), pCross(:,2), pib(:,1), pib(:,2));
                pIn(i, :)   = pCross(in, :);
                pOut(i, :)  = pCross(~in, :);
            end
            do = false;
        catch
            do = true;
        end
    end
end

% -------------------------------------------------------------------------
function [p, t] = addEmpCells(p, t, bnW)
% Add empty cells which appear during the generation of Voronoi grid
    t = sortPtsCounterClockWise(p, t);
    G = tessellationGrid(p, t);
    [fn, pos] = gridFaceNodes(G, (1:G.faces.num));
    fn = reshape(fn, 2, [])';
    assert(all(diff(pos)==2))
    assert(all( fn(:,2) > fn(:,1) ))   
    bnW = [bnW; bnW(1)];
    fW  = zeros(size(bnW,1)-1,1);
    for i = 1 : length(bnW)-1
        n = bnW(i:i+1)';
        n = sort(n);
        fW(i) = find( all(bsxfun(@eq, fn, n), 2) );
    end
    bf  = find( ~all(G.faces.neighbors, 2) );
    bf  = bf(~ismember(bf, fW));
    bfn = fn(bf, :);
    
    for i = 1 : length(fW)
        f = fW(i);
        n = bnW(i:i+1)';
        n = sort(n);
        n1 = n(1);
        n3 = n(2);
        if ~all(G.faces.neighbors(f,:), 2)
           f2 = bf( any(bfn==n1 ,2) );
           [FM, NM] = deal(cell(length(f2),1)); 
           for j = 1 : length(f2)
               n2 = fn(f2(j), :);
               n2 = n2(n2~=n1);
               NM{j} = n2;
               FM{j} = f2(j);
           end
            
           for j = 1 : length(f2)
               while true
                   fm = bf( any(bfn == NM{j}(end), 2) );
                   fm = fm( fm~=FM{j}(end) );
                   nm = fn(fm, :);
                   nm = nm( nm~=NM{j}(end) );
                   NM{j} = [NM{j}, nm];
                   FM{j} = [FM{j}, fm];
                   if any(bnW==nm)
                       break
                   end
               end
           end
           idx = cellfun(@(x)x(end)==n3, NM);
           NM = NM{idx};
           t = [t; {[n1, NM]}];
        end
    end
end

% -------------------------------------------------------------------------
function throwError(L)
    error(['Cannot generate appropriate Voronoi sites, please \n',...
        '   (1) Increase the resolution of well trajectory (add more well points) \n', ...
        'Or (2) Use the grid type ''triangular'' instead \n', ...
        'Or (3) Increase the value of ''WR.ly'', ',...
        'the suggested value is %.0f (may require to enlarge the VOI boundary)\n'], 1.1*max(L));
end
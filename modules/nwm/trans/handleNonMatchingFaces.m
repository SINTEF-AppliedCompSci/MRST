function relation = handleNonMatchingFaces(G1, faces1, G2, faces2, varargin)
% Compute the intersection relation of a surface shared by G1 and G2. 
% The surface will be divided into a set of subfaces due to different face
% gemotries of G1 and G2.
%
% SYNOPSIS:
%   relation = handleNonMatchingFaces(G1, faces1, G2, faces2)
%   relation = handleNonMatchingFaces(G1, faces1, G2, faces2, 'isfaceNodesSorted', true)
%
% PARAMETERS:
%  G1,G2            -  Grids sharing a surface
%  faces1, faces2   -  Surface faces set, faces1 from G1 and faces2
%                      from G2. Basically, faces1 and faces2 constitute the
%                      same 3D surface. The surface is continuous and
%                      completed.
%
% % KEYWORD ARGUMENTS:
%  'isfaceNodesSorted'  -  Wether the nodes of faces stored at 
%                          'G.faces.nodes' are sorted (for both G1 and G2)
%
% RETURNS:
%  relation - Face intersection relation, n x 9 matrix
%              column 1     - Face of G1
%              column 2     - Face of G2
%              column 3     - Areas of intersection subfaces
%              column 4-6   - Centroids of the subfaces
%              column 7-9   - Normals of the subfaces
%
% EXAMPLE:
%     G1 = cartGrid([10, 10, 5], [10, 10, 5]);
%     G2 = cartGrid([20, 20, 5], [10, 10, 5]);
%     G2.nodes.coords = twister(G2.nodes.coords);
%     G2.nodes.coords(:,3) = G2.nodes.coords(:,3) + 5;
%     G1 = computeGeometry(G1);
%     G2 = computeGeometry(G2);
%     figure, hold on
%     plotGrid(G1, 'facecolor', 'none', 'edgecolor', 'r')
%     plotGrid(G2, 'facecolor', 'none', 'edgecolor', 'b'), view(3)
%     faces1 = find(G1.faces.centroids(:,3) == 5);
%     faces2 = find(G2.faces.centroids(:,3) <= 5+0.01 & G2.faces.centroids(:,3) >= 5-0.01);
%     relation = handleNonMatchingFaces(G1, faces1, G2, faces2);
%
% SEE ALSO:
%  `assembleGrids`, `handleMatchingFaces`

    opt = struct('isfaceNodesSorted', false);
    opt = merge_options(opt, varargin{:});
   
    % Nodes of faces
    faceNodes1 = arrayfun(@(f)getSortedFaceNodes(G1, f, opt.isfaceNodesSorted), ...
        faces1, 'UniformOutput', false);
    faceNodes2 = arrayfun(@(f)getSortedFaceNodes(G2, f, opt.isfaceNodesSorted), ...
        faces2, 'UniformOutput', false);

    % Get the intersection relation
    relation = arrayfun(@(f, fn)IntxnRelationSingleFace...
        (G1, f, fn, G2, faces2, faceNodes2), faces1, faceNodes1, ...
        'UniformOutput', false);
    relation = cell2mat(relation);
    
    idx = ~ismember(faces2, relation(:, 2));
    if nnz(idx) > 0
        relation_r = arrayfun(@(f, fn)IntxnRelationSingleFace...
            (G2, f, fn, G1, faces1, faceNodes1), faces2(idx), faceNodes2(idx), ...
            'UniformOutput', false);
        relation_r = cell2mat(relation_r);
        f1 = relation_r(:,2);
        f2 = relation_r(:,1);
        relation_r(:,1) = f1;
        relation_r(:,2) = f2;
        relation = [relation; relation_r];
    end
end

function relation_f = IntxnRelationSingleFace(G1, f1, fn1, G2, faces2, fnodes2)
    % 1. Coordinate transformation ---------------
    fn1 = fn1{1};
    nor_z = G1.faces.normals(f1, :)/G1.faces.areas(f1);
    pts1  = G1.nodes.coords;
    pts2  = G2.nodes.coords;
    [pts1, pts2, T, R] = convertToXYPlane(pts1, fn1, pts2, 'normalZ', []);
    % Closed polygon of face f
    p1 = pts1([fn1; fn1(1)], [1,2]);
    
    check = false;
    if check
        figure(1234), hold on, axis equal off
        coln = 1;
    end
    
    % 2. Find f2 fully located inside f1 ---------
    allIn = cellfun(@(x)all(inpolygon(pts2(x,1), pts2(x,2), p1(:,1), p1(:,2))), fnodes2);
    if nnz(allIn) > 0
        f2   = faces2(allIn);
        R_In = [f2, G2.faces.areas(f2), G2.faces.centroids(f2, :), ...
            G2.faces.normals(f2, :)];
        if check
            plotFaces(G2, f2, 'facecolor', rand(3,1));
            pmid = G2.faces.centroids(f2, :);
            plot3(pmid(:,1), pmid(:,2), pmid(:,3), 's', 'markerfacecolor', 'k')
            drawnow
        end
    else
        R_In = zeros(0, 8);
    end
    
    % 3. Find faces2 intersecting with f1 ---------
    pEdge = addEdgePoints(p1);
    iX = cellfun(@(x)any(inpolygon(pEdge(:,1), pEdge(:,2), ...
        pts2(x,1), pts2(x,2))), fnodes2);
    iX = find(iX);
    if isempty(iX)
        relation_f = [f1*ones(size(R_In,1),1), R_In];
        return
    end
    [f2, areas, f_nor, f_ctd] = deal(cell(length(iX), 1));
    for k = 1 : length(iX)
        fn2   = fnodes2{iX(k)};
        p2    = pts2([fn2; fn2(1)], [1,2]);
        In2   = inpolygon(p1(:,1), p1(:,2), p2(:,1), p2(:,2));
        p1In2 = p1(In2, :);
        In1   = inpolygon(p2(:,1), p2(:,2), p1(:,1), p1(:,2));
        p2In1 = p2(In1, :);
        try
            % Require 'map' toolbox
            [xi, yi] = polyxpoly(p1(:,1), p1(:,2), p2(:,1), p2(:,2), 'unique');
        catch
            % An alternative by 'Gerben J. de Boer'
            [xi, yi] = polyintersect(p1(:,1), p1(:,2), p2(:,1), p2(:,2), 2);
            xyi = unique([xi, yi], 'rows');
            xi = xyi(:,1);
            yi = xyi(:,2);
        end
        p = [p1In2; p2In1; [xi, yi]];
        p = unique(p, 'rows');
        % Remove subfaces whose areas are small compared with areas of both
        % f1 and f2
        minRatio = 0.0001;
        if size(p, 1) > 2
            tmp   = sortPtsCounterClockWise(p, {(1:size(p,1))'});
            p     = p(tmp{1}, :);
            ppoly = [p; p(1,:)];
            area  = polyarea(ppoly(:,1), ppoly(:,2));
            ratio = [area/G1.faces.areas(f1), area/G2.faces.areas(faces2(iX(k)))];
            if all(ratio > minRatio)
                f2{k}    = faces2(iX(k));
                areas{k} = area;
                f_nor{k} = area * nor_z;
                pmid = computeCentroids(p);
                pmid = [pmid, pts1(fn1(1),3)];
                pmid = convertTo3DPlane(pmid, T, R);
                f_ctd{k} = pmid;
                if check
                    ppoly = [ppoly, pts1(fn1(1),3) * ones(size(ppoly,1), 1)];
                    ppoly = convertTo3DPlane(ppoly, T, R);
                    patch( 'Faces',(1:size(ppoly,1)),...
                        'Vertices',ppoly, 'facecolor', rand(3,1));
                    plot3(pmid(:,1), pmid(:,2), pmid(:,3), 's', 'markerfacecolor', 'k')
                    coln = coln + 1;
                    drawnow
                end
            end
        end
    end
    % Remove empty faces
    idx   = ~cellfun(@isempty, areas);
    f2    = f2(idx);
    areas = areas(idx);
    f_ctd = f_ctd(idx);
    f_nor = f_nor(idx);
    R_Ixn = [cell2mat(f2), cell2mat(areas), cell2mat(f_ctd), cell2mat(f_nor)];
    
    % 4. Combine the relations ---------
    relation_f = [R_In; R_Ixn];
    f2 = relation_f(:,1);
    [~, ia] = unique(f2);
    relation_f = relation_f(ia, :);
    relation_f = [f1*ones(size(relation_f,1),1), relation_f];
end


function fn = getSortedFaceNodes(G, f, isSorted)
    fn = G.faces.nodes(G.faces.nodePos(f):G.faces.nodePos(f+1)-1);
    if ~isSorted
        pts1 = G.nodes.coords;
        tmp = zeros(1,3);
        norZ = G.faces.normals(f, :)/G.faces.areas(f);
        pts1 = convertToXYPlane(pts1, fn, tmp, 'normalZ', norZ);
        fn = sortPtsCounterClockWise(pts1(:, [1, 2]), {fn});
        fn = fn{1};
    end
end

function pEdge = addEdgePoints(p)
% Add points along edges of the polygon specified by p
    N = 20;
    pEdge = cell(size(p,1)-1, 1);
    for i = 1 : size(p,1)-1
        p1 = p(i,   :);
        p2 = p(i+1, :);
        d = linspace(0, 1, N+2)';
        pEdge{i} = [ p1(:,1) + d * (p2(:,1)-p1(:,1)), ...
            p1(:,2) + d * (p2(:,2)-p1(:,2)) ];
    end
    pEdge = cell2mat(pEdge);
    pEdge = [p; pEdge];
end

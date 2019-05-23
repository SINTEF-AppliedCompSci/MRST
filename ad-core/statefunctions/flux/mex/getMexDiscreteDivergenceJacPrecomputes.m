function out = getMexDiscreteDivergenceJacPrecomputes(model)
    G = model.G;
    interior = false(G.faces.num, 1);
    interior(model.operators.internalConn) = true;
    facePos = G.cells.facePos;
    cellNo = rldecode(1 : G.cells.num, diff(facePos), 2) .';
    faces = G.cells.faces(:, 1);

    active = interior(faces);

    cellNo = cellNo(active);
    faces = faces(active);

    internal = model.operators.internalConn;
    remap = zeros(G.faces.num, 1);
    remap(internal) = (1:sum(internal))';
    faces = remap(faces);

    N = model.operators.N;
    lc = N(faces, 1);
    rc = N(faces, 2);
    left = (lc ~= cellNo);
    cells = left.*lc + ~left.*rc;
    cells = cells - 1;
    cells(~left) = -cells(~left);



    [c, ps] = rlencode(cellNo);
    if numel(ps) ~= G.cells.num
        % we have some neighborless cells, insert zeros in ps accordingly
        cix = false(G.cells.num, 1);
        cix(cellNo) = true;
        tmp = ps;
        ps = zeros(G.cells.num,1);
        ps(cix) = tmp;
    end
        
    facePos = cumsum([1; ps]);

    localCellIndex = zeros(G.cells.num, 1);
    for i = 1:G.cells.num
        loc = facePos(i):facePos(i+1)-1;
        ca = cells(loc);
        [~, ix] = sort(abs(ca));

        tmp = faces(loc);
        faces(loc) = tmp(ix);
        cells(loc) = ca(ix);

        tmp = find(abs(ca(ix))+1 > i, 1, 'first')-1;
        if isempty(tmp)
            tmp = numel(ca);
        end
        localCellIndex(i) = tmp;

    end
    % Need to sort cells as well
    fpos0 = facePos-1;
    faces = faces - 1;
    out = struct('facePos', fpos0, 'faces', faces, 'cells', cells, 'cellIndex', localCellIndex);
end
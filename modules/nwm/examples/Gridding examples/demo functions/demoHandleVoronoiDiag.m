function [p, t, bdyID] = demoHandleVoronoiDiag(pVor, tVor, pib, pob, pW, tW, bnW)
% An excerpt from 'generateVOIGridNodes' to show how to handle the infinite
% Voronoi diagram. This includes (1) clip the Voronoi diagram, (2) remove
% the conflicting points, and (3) connect with the WR grid.

    % Clip the infinite Voronoi diagram: remove points outside the region
    fdI = @(p)dpoly(p, [pib; pib(1,:)]);
    fdO = @(p)dpoly(p, [pob; pob(1,:)]);
    fd  = @(p)ddiff(fdO(p), fdI(p));
    tol1 = 0.1;
    in  = find(fd(pVor) < tol1 & all(~isinf(pVor), 2));
    map = [(1:length(in))', in];
    pVor = pVor(map(:,2), :);
    
    % Remove conflict points (point too close to each other)
    D = euclideanDistance(pVor, pVor);
    D = triu(D);
    tol2 = 0.1;
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
    
    % Connect with the WR grid: add WR points, and map the connectivity 
    % list again
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


function [c_self, c_other] = boundaryNeighbors(g, subs)
    fa = boundaryFaces(g, subs);
    
    % Extract neighborship
    N = g.faces.neighbors(fa, :);
    % Remove faces on boundary
    bad = any(N == 0, 2);
    fa = fa(~bad);
    N = N(~bad, :);
    
    ok = false(g.cells.num, 1);
    ok(subs) = true;
    ok = [false; ok];
    
    isInt = ok(N(:, 1) + 1);

    [c_self, c_other] = deal(zeros(numel(fa), 1));
    
    c_self( isInt) = N( isInt, 1);
    c_self(~isInt) = N(~isInt, 2);
    
    c_other( isInt) = N( isInt, 2);
    c_other(~isInt) = N(~isInt, 1);

end
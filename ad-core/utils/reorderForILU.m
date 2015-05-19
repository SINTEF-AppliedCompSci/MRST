function [A, b] = reorderForILU(A, b, nc)
    if nargin == 2
        % Include all equations for consideration
        nc = 0;
    end
    d = diag(A);
    d_sub = d((nc+1):end);
    % Find bad entries
    bad = find(d_sub == 0);
    if isempty(bad)
        return
    end
    
    N = size(A, 1);
    
    A_sub = A((nc+1):end, (nc+1):end);
    
    [ii, jj, vv] = find(A_sub);
    
    nbad = numel(bad);
    newInd = zeros(nbad, 1);
    for i = 1:nbad
        % Switch equations around to avoid zero diagonal
        possible = jj(ii == bad(i));
        candidates = ii(jj == bad(i));
        
        new = intersect(possible, candidates);
        if isempty(new)
            warning('Unable to reorder equations, zeros on the diagonal still present');
            continue
        end
        newInd(i) = new(1);
        jj(ii == new(1)) = 0;
    end
    
    renum = (1:N);
    renum(nc + newInd) = nc + bad;
    renum(nc + bad) = nc + newInd;
    
    A = A(renum, :);
    b = b(renum, :);
end

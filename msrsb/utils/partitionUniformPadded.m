function p = partitionUniformPadded(G, dims)
    ijk = gridLogicalIndices(G);
    
    ijk = [ ijk{:} ];
    M = max(ijk) - min(ijk) + 1;

    coarseDim = dims - 1;
    
    
    spacings = cell(G.griddim, 1);
    
    for d = 1:G.griddim
       v = ijk(:, d);
        
       B = coarseDim(d);
       if B == 0
           spacings{d} = M(d);
       else
           cts = lbLinDist((min(v) - 1):max(v) - 1, M(d), B)';
           [~, s] = rlencode(cts);

           divEl = s(1);
           e1 = floor(divEl/2);
           e2 = divEl - e1;
           spacings{d} = [e1; s(2:end); e2];
       end
    end
    p = partitionTensor(G, spacings{:});
end

function f = lbLinDist(f, M, B)
    % Stolen from partitionUI

    L = floor(M ./ B);  % Tentative number of cells per coarse block.
    R = mod(M, B);      % Additional cells not previously accounted for.
    f = max(floor(f ./ (L + 1)), floor((f - R) ./ L));
end
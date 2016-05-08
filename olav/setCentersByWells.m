function CG = setCentersByWells(CG, W)
    for i = 1:numel(W)
        c = W(i).cells;
        local = unique(CG.partition(c));
        for j = 1:numel(local)
            cix = c(CG.partition(c) == local(j));
            mid = ceil(numel(cix)/2);
            CG.cells.centers(local(j)) = cix(mid);
        end
    end
end

function CGm = changeCoarseCenterToWellLocation(CGm,W)
% Designates coarse nodes to cells containing wells. For coarse blocks with
% multiple cells, it selects a random well location as coarse node
if CGm.parent.griddim>2
    CGm = setCentersByWells(CGm, W);
else
    wcells = [W.cells];
    for i = 1:numel(wcells)
        in = CGm.partition(wcells(i));
        CGm.cells.centers(in) = wcells(i);
    end
end
return

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
return
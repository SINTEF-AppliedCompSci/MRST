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
    
function CG = extendInteractionToFineFaces(CG, coarsefaces)
%     coarseFaceNo = rldecode((1:CG.faces.num)', diff(CG.faces.connPos));
    
%     coarsefaces = unique(coarseFaceNo(finefaces));
    coarsefaces = unique(coarsefaces);
    
    isBF = false(CG.faces.num, 1);
    isBF(boundaryFaces(CG)) = true;
    
    coarsefaces = coarsefaces(isBF(coarsefaces));
    
    facecells = sum(CG.faces.neighbors(coarsefaces, :), 2);
    
    for i = 1:numel(coarsefaces)
        cf = coarsefaces(i);
        
        ci = facecells(i);
        
%         n = neighbors(CG, facecells(i));
        candidates = coarseNeighbors(CG, ci, false);
        candidates = setdiff(candidates, facecells);
        for j = 1:numel(candidates)
            CG.cells.interaction{candidates(j)} =...
                vertcat(CG.cells.interaction{candidates(j)}, ...
                        CG.cells.interaction{ci});

%             CG.cells.interaction{candidates(j)} =...
%                 vertcat(CG.cells.interaction{candidates(j)}, find(CG.partition == facecells(i)));
        end
    end
    
    for i = 1:CG.cells.num
        CG.cells.interaction{i} = unique(CG.cells.interaction{i});
    end
%     indicator = fal
end

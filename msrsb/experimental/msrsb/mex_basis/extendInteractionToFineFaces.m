function CG = extendInteractionToFineFaces(CG, coarsefaces)
%Undocumented Utility Function

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

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

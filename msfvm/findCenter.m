function [index localindex] = findCenter(cg, center, blockInd)
%indices of blocks in the current cell
%blockInd = find(cg.partition == coarseindex);
%find the centroids

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


    blockCentroids = cg.parent.cells.centroids(blockInd,:);
    if cg.griddim == 3
        %find distance from coarse centroid to fine centroids
        centroidDist =  (blockCentroids(:,1)-center(1)).^2 + ...
                        (blockCentroids(:,2)-center(2)).^2 + ...
                        (blockCentroids(:,3)-center(3)).^2;
    else
        centroidDist =  (blockCentroids(:,1)-center(1)).^2 + ...
                        (blockCentroids(:,2)-center(2)).^2;
    end
    %find the index of the least distance
%     least = find(centroidDist == min(centroidDist),1);
    [val least] = min(centroidDist); %#ok MATLAB requires this output
    %lookup global index and take "first" element in case several match
    index = blockInd(least(1));
    localindex = least(1);
end

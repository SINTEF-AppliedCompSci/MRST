function [xMin, xMax, ixMin, ixMax, neg, pos] = getMinMax(x, nn)
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

    % Prepare sparse matrix vectors
    num = numel(nn);
    ii  = [rldecode((1:num)', nn, 1);
           rldecode((1:num)', max(nn) - nn, 1)];
    jj  = [mcolon(1, nn)';
           mcolon(nn+1, max(nn))'];
    nanVec = nan(sum(max(nn) - nn), 1);
    
    % Compute minimum and maximum cell coordinates
    [xMin, xMax, ixMax, ixMin, neg, pos] = deal(zeros(num, size(x,2)));
    for dNo = 1:size(x,2)
        v = [x(:,dNo); nanVec];
        xMat = sparse(ii, jj, v, num, max(nn));
        [xMin(:,dNo), ixMin(:,dNo)] = min(xMat, [], 2);
        [xMax(:,dNo), ixMax(:,dNo)] = max(xMat, [], 2);
        neg(:,dNo) = all(xMat < 0 | isnan(xMat),2);
        pos(:,dNo) = all(xMat > 0 | isnan(xMat),2);
    end
    
end

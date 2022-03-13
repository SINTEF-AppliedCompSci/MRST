function [p, binAngle, binRadius, dist] = refineNearWell(pts, wellpt, varargin)
%Partition a set of points based on proximity to a well point
%
% SYNOPSIS:
%   p = refineNearWell(points, wellpoint)
%
% REQUIRED PARAMETERS:
%   pts    - Set of points that are to be partitioned (for example
%            G.cells.centroids)
%
%   wellpt - A single point to be used for partitioning (for example the
%            well heel. Routine currently only partitions in xy-plane.
%
% OPTIONAL PARAMETERS:
%  maxRadius   - The distance in meters for which a partition will be
%                applied. (Defaults to inf, i.e. all points). If maxRadius
%                is a two-vector, the partition is applied inside a
%                rectangle rather than inside a circle.
%
%  angleBins   - The number of sectors the region will be divided into
%                along the angle around the point. 
%
%  radiusBins  - Number of distinct regions produced as the radius
%                increases.
%  
%  logBins     - If true (default) then the radius bins will have width
%                that drops off as a log(r). Otherwise, will distribute in
%                a linear fashion.
%
% RETURNS:
%   p          - Partition vector. One entry for each point. Will contain
%                zero for any points that are out of bounds for the
%                maxRadius parameter.
%
%  binAngle   - The partition used for the sector partitioning.
%
%  binRadius  - The partition used for the radius partitioning.
%
%  dist       - Distance per point used for cutoff by maxRadius.
%
%  
%
% SEE ALSO:
% `partitionUI`
%   

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


    opt = struct('maxRadius',      inf,  ...
                 'angleBins',      6, ...
                 'radiusBins',     3, ...
                 'logBins',        true ...
                 );
    
    opt = merge_options(opt, varargin{:});
    
    assert(numel(opt.radiusBins) == 1);
    assert(numel(opt.angleBins) == opt.radiusBins || numel(opt.angleBins) == 1, ...
        ['Vector of angular bins must either be either one value for each radial', ...
        ' bin, or one value that will be used for all']);
    
    nbins = numel(opt.angleBins);
    npts  = size(pts, 1);
    
    % Paritition the angles (once for each radial refinement index)
    binAngle = zeros(npts, nbins);
    for i = 1:numel(opt.angleBins)
        binAngle(:, i) = partitionSectors(pts, wellpt, opt.angleBins(i));
    end
    % Partition into bins along the increasing radius
    [binRadius, dist, outside] = partitionBins(pts, wellpt, ...
                                opt.radiusBins, opt.logBins, opt.maxRadius);
    
    
    % Use both sector and distance indicators
    if nbins == 1
        p = binRadius + max(binRadius).*binAngle;
    else
        p = nan(npts, 1);
        for i = 1:nbins
            loc = binRadius == i;
            p(loc) = binRadius(loc) + max(binRadius).*binAngle(loc, i);
        end
    end
    % Ensure partition does not skip numbers
    p = compressPartition(p);
    % Honor max distance from center. Mask to zero, so partitions can be
    % created simply by addition.
    p(outside) = 0;
end


function [ind, dist, outside] = partitionBins(pts, p0, binNum, useLogBins, maxRadius)
    if numel(maxRadius)==1
       dist = len(bsxfun(@minus, pts(:, 1:2), p0(1:2)));
       outside = dist > maxRadius;
    else
       dist = abs(bsxfun(@minus, pts(:,1:2), p0(1:2)));
       outside = any(bsxfun(@rdivide, dist, maxRadius)>1,2);
       dist = len(dist);
   end
    d = normalize(dist);
    if useLogBins
        d = d*10 + 1;
        d = log10(d);
    end
    d(outside) = max(d(~outside));
    ind = groupByBins(d, binNum);
end

function ind = partitionSectors(pts, p0, sectornum)
    phi = getAngleTan(p0, pts);
    ind = groupByBins(phi, sectornum, -pi, pi);
end

function phi = getAngleTan(p0, p1)
    p0 = repmat(p0, size(p1, 1)/size(p0, 1), 1);
    
    v = p1 - p0;
    v = bsxfun(@rdivide, v, len(v));

    phi = atan2(v(:, 1), v(:, 2));
end

function l = len(v)
    l = sqrt(sum(v.^2, 2));
end

function ind = groupByBins(data, numBins, dMin, dMax)
    if numBins == 1
        ind = ones(size(data));
        return
    end
    if nargin < 3
        dMin = min(data);
    end
    if nargin < 4
        dMax = max(data);
    end
    data = data - dMin;
    data = data./(dMax - dMin);
    
    ind = ones(size(data, 1), 1);
    w = 1/numBins;
    for i = 1:numBins
        int = data > (i-1)*w & data <= i*w;
        ind(int) = i;
    end
end

function data = normalize(data)
    data = data - min(data);
    data = data./max(data);
end

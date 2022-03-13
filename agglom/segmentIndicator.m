function [partition, edges] = segmentIndicator(G, indicator, bins, varargin)
%Segments a fine grid into blocks according to indicator.
%
% SYNOPSIS:
%   partition = segmentIndicator(G, indicator, bins)
%   partition = segmentIndicator(G, indicator, bins, 'pn1', pv1, ...)
%   [partition, edges] = segmentIndicator(...)
%
% DESCRIPTION:
%   This function segments fine grid cells into prescribed bins according
%   the given indicator value (indicator).  The cell groupings may be split
%   into connected components.
%
% REQUIRED PARAMETERS:
%   G         - Grid data structure discretising the reservoir model
%               (fine grid, geological model).
%
%   indicator - Cell-wise value of some measure/indicator function.
%               Assuming positive numbers.
%
%   bins      - Gives the bins to segment the fine grid cells into if it is
%               a vector and the number of bins if a scalar.
%
% OPTIONAL PARAMETERS:
%
%   split    - Whether or not to split the grouped cells into connected
%              components. Default: true
%
%   verbose  - Whether or not display number of blocks in the resulting
%              partition. Default value dependent upon global verbose
%              settings of function 'mrstVerbose'.
%
%
% RETURNS:
%   partition - Partition vector after segmenting all cells into blocks
%               according to the indicator
%
%   edges     - Segmentation levels.  Equal to 'bins' if a vector.
%               otherwise, in case of scalar, positive 'bins',
%
%               edges = LINSPACE(MIN(indicator), MAX(indicator), bins + 1)
%
% SEE ALSO:
%   `mergeBlocks`, `refineBlocks`

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


opt  = struct('verbose', mrstVerbose, 'split', true);
opt = merge_options(opt, varargin{:});

if numel(bins)==1
   Imin = min(indicator);
   Imax = max(indicator);
   I = (indicator - Imin)/(Imax - Imin + eps);
   partition = round( I*bins + 0.5);
   edges = linspace(Imin, Imax, bins + 1);
else
   % Use the built-in function histc to perform the segmentation. Infinity
   % is added as the last edge value for the bins to avoid the maximum
   % value of the indicator to end up in its own bin.
   edges = bins;
   edges(end) = inf;
   edges = unique(edges);

   [N, partition] = histc(indicator, edges); %#ok<ASGLU>
end

if opt.split
   require coarsegrid
   partition = compressPartition(processPartition(G, partition));
end
dispif(opt.verbose, 'SegmentIndicator: %d blocks\n', max(partition));

end

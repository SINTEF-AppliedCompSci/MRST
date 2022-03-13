function p = mexPartitionMETIS(G, n, varargin)
%Partition Graph/Grid Using METIS Library
%
% SYNOPSIS:
%   p = mexPartitionMETIS(G, numBlocks)
%   p = mexPartitionMETIS(G, numBlocks, w)
%   p = mexPartitionMETIS(G, numBlocks, w, 'pn1', pv1, ...)
%
% PARAMETERS:
%   G         - Grid structure.
%
%   numBlocks - Number of blocks into which to partition the grid 'G'.
%
%   w         - Connection weighting array.  Size equal to number of
%               half-faces (size(G.cells.faces, 1), number of faces
%               (G.faces.num), or number of internal faces
%               (sum(all(G.faces.neighbors > 0, 2)).  Often just the
%               one-sided transmissibilities calculated by function
%               `computeTrans` or the total face transmissibilities
%               calculated by function `getFaceTransmissibility`.
%
% OPTIONAL PARAMETERS:
%   ufactor - Percentage overweight allowed in coarse blocks. Examples:
%                   1.01: Largest difference in cells per coarse block is
%                         equal to or less than 1%.
%
%                   1.50: Largest difference in cells per coarse block is
%                         equal to or less than 50%. (DEFAULT)
%
%   no2hop  - Do not use multiple hops when coarsening. Passed on to METIS.
%
%   seed    - METIS seeding parameter to create consistent results for
%             given platform.  Default: 0.
%
%   ncuts   - Passed on to METIS.
%
%   useLog  - Log_10 transform weighting array `w` before creating
%             connection strengths.  This can be useful if the weighting is
%             based on a transmissibility field which varies by several
%             orders of magnitude locally.
%
% RETURNS:
%   p - Partition vector.  Contains no empty or multiply connected blocks.
%
% COMPATIBILITY:
%   Tested on MATLAB R2020a
%     - Windows 10, 1909, Visual Studio 15.9.27 (MSVC 19.16.27043)
%     - Ubuntu Linux 18.04.5 LTS, GCC 7.5.0
%
% SEE ALSO:
%   `partitionMetis`, `compressPartition`, `processPartition`.

%{
Copyright 2020-2022 SINTEF Digital, Mathematics & Cybernetics.

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

   opt = struct('ufactor',     1.5, ...
                'ncuts',       1,   ...
                'seed',        0,   ...
                'niter',      10,   ...
                'minconn',  true,   ...
                'useLog',  false,   ...
                'no2hop',  false);

   w = [];
   if (nargin > 2) && (mod(nargin, 2) == 1)
      % p = mexPartitionMETIS(G, n, w, ...)
      w = varargin{1};
      varargin = varargin(2 : end);
   end

   opt = merge_options(opt, varargin{:});

   if ~isempty(w) && opt.useLog
      w = log10(w);
   end

   mexFuncArgs = { opt };
   if ~isempty(w)
      mexFuncArgs = [ mexFuncArgs, { w } ];
   end

   p = mexPartitionMETIS_Impl(G, n, mexFuncArgs{:});
end

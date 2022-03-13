function [p, A] = partitionMETIS(G, T, n, varargin)
%Partition grid according to connection strengths
%
% SYNOPSIS:
%    p     = partitionMETIS(G, T, n)
%    p     = partitionMETIS(G, T, n, 'pn1', pv1, ...)
%   [p, A] = partitionMETIS(...)
%
% PARAMETERS:
%   G - Grid structure.
%
%   T - One-sided transmissibilities as defined by function 'computeTrans'.
%       Possibly affected by mobility.
%
%   n - Number of blocks into which to partition the grid 'G'.
%
% OPTIONAL PARAMETERS:
%   ufactor - Percentage overweight allowed in coarse blocks. Examples:
%                   1.01: Largest difference in cells per coarse block is
%                         equal to or less than 1%.
%
%                   1.50: Largest difference in cells per coarse block is
%                         equal to or less than 50%. (DEFAULT)
%
%   no2hop  - Do not use multiple hops when coarsening. Passed onto METIS.
%
%   seed    - METIS seeding parameter to create consistent results for
%             given platform. Default: 0.
%
%   ncuts   - Passed onto METIS.
%
%   useLog  - Log_10 transform transmissibilities before creating
%             connection matrix. This can be useful when transmissibilities
%             vary by several orders of magnitude locally.
%
% RETURNS:
%   p - Partition vector.  Contains no empty or multiply connected blocks.
%
%   A - Matrix used for partitioning.
%
% SEE ALSO:
%   `incompTPFA`, `callMetisMatrix`, `compressPartition`, `processPartition`.

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

   require coarsegrid incomp

   opt = struct('ufactor', 1.5, ...
                'ncuts',   1, ...
                'seed',    0, ...
                'niter',  10, ...
                'minconn', true, ...
                'useLog',  false, ...
                'no2hop',  false);
   opt = merge_options(opt, varargin{:});

   fluid = initSingleFluid('mu', 1, 'rho', 0*kilogram/meter^3);

   nr = ceil(n);
   if nr ~= n
       warning(['Number of blocks should be a natural number! Rounding up: ', ...
           num2str(n), ' -> ', num2str(nr), '.']);
      n = nr;
   end
   % Assemble matrix of connection strengths.
   %
   % Note: Use a no-op linear solver to avoid costly linear solve that's
   % not actually needed in this case.
   if opt.useLog
       T = log10(T);
   end
   x = incompTPFA(initState(G, [], 0), G, T, fluid, ...
                  'LinSolve', @(A, x) zeros(size(x)), ...
                  'use_trans',  numel(T) == G.faces.num, ...
                  'MatrixOutput', true);
   A = x.A;
   % Call METIS with the
   %  - 'minconn' option (favours relatively square, convex blocks in
   %    regions of low contrast variance)
   %  - 'contig' to ensure contiguous coarser grids
   %  - 'ufactor' set to 500 to allow for 150% variation on coarse block
   %    size
   %  - 'cut' is default objtype
   opts = sprintf('-contig -ufactor=%d -objtype=cut -ncuts=%d -seed=%d -niter=%d', ...
                  floor(1000 * (opt.ufactor - 1)), opt.ncuts, opt.seed, opt.niter);
   if opt.no2hop
       opts = [opts, ' -no2hop'];
   end
   if opt.minconn
       opts = [opts, ' -minconn'];
   end
   p = callMetisMatrix(A, n, opts);
   p = processPartition(G, compressPartition(p));
end

function src_cg = convertSource2Coarse(CG, src)
%Accumulate fine-scale source terms to coarse scale
%
% SYNOPSIS:
%   src_cg = convertSRC2Coarse(CG, src)
%
% PARAMETERS:
%   CG  - Coarse grid as defined by function 'generateCoarseGrid'.
%
%   src - Source structure as defined by function 'addSource'.
%
% RETURNS:
%   src_cg - Source data structure defined on coarse grid.
%
% SEE ALSO:
%   `convertBC2Coarse`, `addSource`.

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


   if isempty(src) src_cg = []; return; end

   pq          = CG.partition(src.cell);
   has_src     = false([CG.cells.num, 1]);
   has_src(pq) = true;

   has_pos = accumarray(pq, src.rate > 0, [CG.cells.num, 1]);
   has_neg = accumarray(pq, src.rate < 0, [CG.cells.num, 1]);

   assert (~any(has_pos & has_neg), ...
          ['Unsafe configuration: Injection and production in ', ...
           'same coarse block']);

   blk = find(has_src);

   blk_rate = accumarray(pq, src.rate, [CG.cells.num, 1]);

   if isfield(src, 'sat'),
      sat = zeros([CG.cells.num, size(src.sat, 2)]);
      for i = 1 : size(sat,2),
         sat(:,i) = accumarray(pq, src.rate .* src.sat(:,i), ...
                               [CG.cells.num, 1]) ./ blk_rate;
      end

      sat_args = { 'sat', sat(blk, :) };
   else
      sat_args = {};
   end

   src_cg   = addSource([], blk, blk_rate(blk), sat_args{:});
end

function Wc = coarsenWells(Gc, W)
% Generate well structure according to coarse grid Gc
% SYNOPSIS:
%   Wc = coarsenWells(Gc, W)
%
% PARAMETERS:
%   Gc       - coarse grid 
%   W        - Well structure compatible with Gc.parent 

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

   Wc = W;
   for i = 1 : numel(W)
      fcells  = W(i).cells;
      cgcells = Gc.partition(fcells);

      tab        = sortrows([cgcells, fcells, (1 : numel(fcells)) .']);
      [cells, n] = rlencode(tab(:,1));
      fcellspos  = cumsum([1 ; n]);

      if Gc.griddim > 2
         pno = rldecode(1 : numel(cells), n, 2) .';
         cc  = Gc.parent.cells.centroids(tab(:,2), 3);
         cv  = Gc.parent.cells.volumes  (tab(:,2));

         hpos = sparse(pno, 1 : numel(pno), cv) ...
                * [ cc, ones([numel(pno), 1]) ];

         hpos = hpos(:,1) ./ hpos(:,2);         clear pno cc cv
      else
         hpos = 0;
      end

      Wc(i).cells     = cells;
      Wc(i).WI        = nan([numel(cells), 1]);
      Wc(i).dZ        = hpos - W(i).refDepth;
      Wc(i).fcellspos = fcellspos;
      Wc(i).fcells    = tab(:,2);
      Wc(i).fperf     = tab(:,3);
      Wc(i).parent = W(i);
   end
end

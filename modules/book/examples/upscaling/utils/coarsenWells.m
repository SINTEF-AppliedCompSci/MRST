function Wc = coarsenWells(Gc, W)
% Generate well structure according to coarse grid Gc
% SYNOPSIS:
%   Wc = coarsenWells(Gc, W)
%
% PARAMETERS:
%   Gc       - coarse grid
%   W        - Well structure compatible with Gc.parent
   Wc = W;
   for i = 1 : numel(W),
      fcells  = W(i).cells;
      cgcells = Gc.partition(fcells);

      tab        = sortrows([cgcells, fcells, (1 : numel(fcells)) .']);
      [cells, n] = rlencode(tab(:,1));
      fcellspos  = cumsum([1 ; n]);

      if Gc.griddim > 2,
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

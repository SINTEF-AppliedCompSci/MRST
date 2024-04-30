function cgwells  = makeCGwells(cg,wells)
%Undocumented Utility Function

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

    cgwells=wells;
    for i=1:numel(wells)
       fcells= wells(i).cells;
       cgcells=cg.partition(fcells);
       %nperf=cellfun('prodofsize',[wells.cells]);
       %wno=rldecode(1:numel(wells),nperf);
       tab=sortrows([cgcells,fcells,[1:numel(fcells)]']);%#ok
       [cells,n]=rlencode(tab(:,1));
       fcellspos =cumsum([1;n]);
       pno=rldecode([1:numel(cells)]',diff(fcellspos));%#ok
       if(cg.griddim>2)
        hpos=accumarray(pno,cg.parent.cells.centroids(tab(:,2),3)...
              .*cg.parent.cells.volumes(tab(:,2)))./...
            accumarray(pno,cg.parent.cells.volumes(tab(:,2)));
       else
           hpos=0;
       end
       cgwells(i).cells=cells;
       cgwells(i).WI=nan(numel(cells),1);
       cgwells(i).dZ=hpos-wells(i).refDepth;
       cgwells(i).fcells=tab(:,2);
       cgwells(i).fcellspos=fcellspos;
       cgwells(i).fperf=tab(:,3);
       %{
       cgwells=[cgwells,struct('cells'   , cells,           ...
             'type'    , wells(i).type,             ...
             'val'     , wells(i).val,              ...
             'WI'      , nan(numel(cells),1),                   ...
             'dZ'      , hpos-wells(i).refDepth, ...
             'name'    , wells(i).name,             ...
             'compi'   , wells(i).compi,           ...
             'refDepth', wells(i).refDepth,         ...
             'sign'    , wells(i).sign,...
             'fcells'  , tab(:,2),...%not needed???
             'fcellspos',fcellspos,...
             'fperf',tab(:,3))];%#ok
       %}
   end
   cgwells(i).parent=wells(i);
end

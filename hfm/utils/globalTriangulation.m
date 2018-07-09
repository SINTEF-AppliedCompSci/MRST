function GlobTri = globalTriangulation(G)
% This function creates a global triangulation such that each cell in a
% grid is made up of one or more full tetrahedrons. The global
% triangulation can then be used to locate the cells that enclose a given
% set of points.

%{
Copyright 2009-2015: TU Delft and SINTEF ICT, Applied Mathematics.

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

[cn,cnmap] = gridCellNodes(G,1:G.cells.num); i = 1;

tx = delaunayTriangulation(G.nodes.coords(cn(cnmap(i):cnmap(i+1)-1),:));
P = tx.Points;
T = tx.ConnectivityList;
map = [];
map(1:size(T,1),1) = 1;
lookup = zeros(G.nodes.num,1);
lookup(cn(cnmap(i):cnmap(i+1)-1)) = transpose(1:size(P,1));
for i = 2:G.cells.num
    nx = cn(cnmap(i):cnmap(i+1)-1);
    tx = delaunayTriangulation(G.nodes.coords(nx,:));
    map(end+1:end+size(tx,1)) = i;
    if any(lookup(nx))
        conn = tx.ConnectivityList;
        loc = lookup(nx);
        local_ind  = find(loc);
        conn2 = zeros(size(conn));
        for j = 1:numel(local_ind)
            conn2(conn == local_ind(j)) = loc(local_ind(j));
        end
        local_add = find(loc==0);
        global_add = size(P,1) + 1:size(P,1) + numel(local_add);
        for j = 1:numel(local_add)
            conn2(conn == local_add(j)) = global_add(j);
        end
        lookup(nx(local_add)) = global_add;
        P = [P;G.nodes.coords(nx(loc==0),:)]; %#ok
        T = [T;conn2]; %#ok
    else
        P = [P;G.nodes.coords(nx,:)]; %#ok
        conn = tx.ConnectivityList + size(P,1);
        T = [T;conn]; %#ok
    end
end
Gtri = triangulation(T,P);
GlobTri = struct;
GlobTri.Tri = Gtri;
GlobTri.map = map;
return
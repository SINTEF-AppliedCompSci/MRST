function GC = generateHybridCoarseGrid(G, partition)
%Undocumented Utility Function

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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

    GC = generateCoarseGrid(G, partition);
    GC = coarsenGeometry(GC);
    
    faces = GC.cells.faces(:,1);
    faces = GC.faces.fconn(mcolon(GC.faces.connPos(faces), GC.faces.connPos(faces+1)-1));
    
    nff = diff(GC.faces.connPos);
    ncf = accumarray(rldecode((1:GC.cells.num)', diff(GC.cells.facePos), 1), nff(GC.cells.faces(:,1)));
    GC.cells.faces = faces;
    
    GC.cells.facePos = [0; cumsum(ncf)] + 1;
    
    GC.faces = G.faces;
    GC.faces.neighbors = G.faces.neighbors;
    ix = GC.faces.neighbors ~= 0;
    GC.faces.neighbors(ix) = partition(G.faces.neighbors(ix));
    GC.faces.fconn = (1:GC.faces.num)';
    GC.faces.connPos = (1:GC.faces.num+1)';
    
    GC.edges = G.edges;
    
    GC.nodes = G.nodes;
    
    GC.type = [G.type, { mfilename }];
    
end

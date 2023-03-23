function frac = getIndepNetwork(fl)
% getIndepNetwork(fl) extracts independant fracture networks from a set of
% 'n' fracture lines given by an n-by-4 matrix with each row containing the
% end points of a fracture line. This function is used when the matrix is
% two-dimensional.
%
% SYNOPSIS:
%   frac = getIndepNetwork(fl)
%
% REQUIRED PARAMETERS:
%   fl - fracture lines represented by it's end points as [x1 y1 x2 y2]. fl
%        will have 1 row per fracture line.
%
% RETURNS:
%   frac - Structure containing information about individual fractures
%          with the following sub-structures: (a) lines - Structure with
%          the fields "network"
%                      (network to each fracture belongs), "endp"
%                      (endpoints of each fracture line as supplied by
%                      'fl') and "cells" (matrix cells containing the
%                      fracture) . Size = 1-by-rows(fl).
%          (b) network - stores indices for fracture lines contained in
%                        each network.
% SEE ALSO:
%   markcells2D

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


out = lineSegmentIntersect(fl,fl);
if ~isequal(out.intAdjacencyMatrix,out.intAdjacencyMatrix.')
    warning('Some intersections not registered. Fracture networks may not be entirely independant');
end
network = struct;
count = 0;
flines = 1:size(fl,1);
while(~isempty(flines))
    count = count+1;
    addf = flines(1);
    network(count).lines = [];
    conn2 = find(out.intAdjacencyMatrix(addf,:));
    conn = addf;
    while(~isempty(conn2))
        conn = [conn,conn2]; %#ok
        conn_temp = [];
        for i = 1:numel(conn2)
            addf = find(out.intAdjacencyMatrix(conn2(i),:));
            addf = addf(~ismember(addf,conn));
            conn_temp = [conn_temp, addf]; %#ok
        end
        conn2 = conn_temp;
    end
    conn = unique(conn);
    network(count).lines(end+1:end+numel(conn)) = conn;
    flines = flines(~ismember(flines,conn));
end
frac.network = network;                
for i = 1:numel(frac.network)
    for j = 1:numel(frac.network(i).lines)
        frac.lines(frac.network(i).lines(j)).network = i;
        frac.lines(frac.network(i).lines(j)).endp = fl(frac.network(i).lines(j),:);
    end
end
return
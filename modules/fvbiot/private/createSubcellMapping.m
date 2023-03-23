function [cno, nno, hfno, fno, subfno, subhfno] = createSubcellMapping(g)
% Create numbering of sub-cell quantities needed for mpsa/mpfa
% discretization.
%
% This method is modified from a file in The MATLAB Reservoir Simulation
% Toolbox (MRST), see the sub-function createMapping within 
%   mrst/modules/mpfa/computeMultiPointTrans.m
% 
%{
Parital copyright 2009-2016 SINTEF ICT, Applied Mathematics.
Partial copyright 2016, University of Bergen.

This file is part of FVBiot.

FVBiot is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

FVBiot is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%}

cellno   = rldecode(1:g.cells.num, diff(g.cells.facePos), 2) .';
col      = 1+(cellno == g.faces.neighbors(g.cells.faces(:,1), 2));
nhfaces  = g.cells.facePos(end)-1;
hfaces   = accumarray([g.cells.faces(:,1), col], 1:nhfaces);
hfaces   = rldecode(hfaces, diff(g.faces.nodePos));


cells    =  rldecode(g.faces.neighbors, diff(g.faces.nodePos));
nodes    =  repmat(g.faces.nodes, [2,1]);
faces    =  repmat(rldecode(1:g.faces.num, diff(g.faces.nodePos),2)', [2,1]);
subfaces =  repmat((1:size(g.faces.nodes,1))', [2,1]);
i        =  cells~=0;
w        =  [cells(i), nodes(i), hfaces(i), faces(i), subfaces(i)];
w        =  double(sortrows(w));


cno     = w(:,1);
nno     = w(:,2);
hfno    = w(:,3);
fno     = w(:,4);
subfno  = w(:,5);
subhfno = (1:numel(cno))';
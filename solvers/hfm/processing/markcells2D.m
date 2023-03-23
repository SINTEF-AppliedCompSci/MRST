function [G, frac] = markcells2D(G,frac)
% markcells2D assigns a indicator to cells containing embedded fractures or
% discrete fractures at any of its faces. Indices of fracture lines and
% networks are also stored for the corresponding cells. Each fracture is
% assumed to be a line by this function. See processFracture2D for details
% on input and output.

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


G.cells.fracture.indicator = zeros(G.cells.num,1);
G.cells.fracture.dfm_indicator = zeros(G.cells.num,1);
G.cells.fracture.dfm_line_num = cell(G.cells.num,1);
G.cells.fracture.line_num = cell(G.cells.num,1);
G.cells.fracture.network = cell(G.cells.num,1);
for j = 1:numel(frac.lines)
    frac.lines(j).cells = [];
end
xy = struct;
csize = double(G.griddim<3)*sqrt(min(G.cells.volumes)) + ...
        double(G.griddim>3)*(min(G.cells.volumes))^(1/3);
for j = 1:numel(frac.lines)
    l1 = frac.lines(j).endp;
    xe = l1(1:2:3); ye = l1(2:2:4);
    r = csize/(sqrt(diff(xe)^2 + diff(ye)^2))/10;
    t = transpose(0:r:1);
    t = [t(1:end-1);unique([t(end);1])];
    xe = xe(1)*t+xe(2)*(1-t);
    ye = ye(1)*t+ye(2)*(1-t);
    xy(j).xe = xe;
    xy(j).ye = ye;
end
%% Check for faces intersecting with fracture line segments
for i = 1:G.cells.num
    faces = G.cells.faces(G.cells.facePos(i):G.cells.facePos(i+1)-1,1);
    cnodes = [];
    for k = 1:numel(faces)
        fnodes = G.faces.nodes(G.faces.nodePos(faces(k)):G.faces.nodePos(faces(k)+1)-1); % 2 for 2D, 4 for 3D
        cnodes = [cnodes;fnodes]; %#ok
    end
    pts = G.nodes.coords(cnodes,:);
    k = convhull(pts,'simplify',true);
    for j = 1:numel(frac.lines)
        xx = xy(j).xe; yy = xy(j).ye;
        [in,on] = inpolygon(xx,yy,pts(k,1),pts(k,2)); 
        if any(in-on)
            frac.lines(j).cells(end+1) = i;
            G.cells.fracture.indicator(i) = 1;
            G.cells.fracture.line_num{i,1} = ...
                [G.cells.fracture.line_num{i,1},j];
            G.cells.fracture.network{i,1} = ...
                [G.cells.fracture.network{i,1},frac.lines(j).network];
        elseif numel(find(on))>2 && all(in==on) 
            frac.lines(j).cells(end+1) = i;
            G.cells.fracture.dfm_indicator(i) = 1; % dfm type fracture on edges
            G.cells.fracture.dfm_line_num{i,1} = ...
                [G.cells.fracture.dfm_line_num{i,1},j];
            G.cells.fracture.network{i,1} = ...
            [G.cells.fracture.network{i,1},frac.lines(j).network];
        end
    end
end

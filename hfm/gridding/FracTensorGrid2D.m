function G = FracTensorGrid2D(G,F,a,varargin)
% FracTensorGrid2D uses information about each fracture line supplied by
% 'F' obtained through assembleFracNodes2D and creates a tensorGrid for
% each fracture line using the fracture aperture 'a'.
%
% SYNOPSIS:
%   G = FracTensorGrid2D(G,F,a)
%   G = FracTensorGrid2D(G,F,a,'plot')
%
% REQUIRED PARAMETERS:
%
%   G  - Matrix grid data structure.
%
%   F  - Structure containing information abo'Frac'ut partitioned fracture
%        lines as returned by assembleFracNodes2D.
%
%   a  - Fracture aperture, Must either have a single value for all
%        fracture lines or an individual value for each fracture line.
%        Varying fracture apartures for one line are not supported.
%
% OPTIONAL PARAMETER:
%
%   'plot' - Plot all fracture grid structures returned under G.FracGrid.
%
% RETURNS:
%   G - Grid structure with sub-structure G.FracGrid. G.FracGrid in-turn
%       contains a grid structure for each fracture line with nomenclature
%       such as Frac1, Frac2, and so on. G.FracGrid.Frac1, for example,
%       will have the same basic grid structure as G.
%
% SEE ALSO:
%   assembleFracNodes2D, tensorGrid

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

if numel(a) == 1, a = repmat(a,numel(F),1); 
else assert(numel(a) == numel(F),'Either specify 1 aperture per fracture line or 1 for all'); end

Gf = struct();
face_start = G.faces.num+1;
for i = 1:numel(F)
    fieldname = ['Frac',num2str(i)];
    coords = F(i).nodes.coords;
    endp = [coords(1,:);coords(end,:)];
    diff1 = diff(endp,1);
    if abs(diff1(2)) < eps*100 % // to x-axis
        temp = tensorGrid(coords(:,1),unique(roundsd([coords(1,2);coords(1,2)+a(i)],14),'rows'));
        Gf.(fieldname) = computeGeometry(temp);
    elseif abs(diff1(1)) < eps*100 % // to y-axis
        temp = tensorGrid(unique(roundsd([coords(:,1); coords(:,1)+a(i)],14),'rows'),coords(:,2));
        Gf.(fieldname) = computeGeometry(temp);
    else
        new_coords = translateLine(coords,a(i));
        m = diff1(2)/diff1(1); % slope
        theta = atand(m);
        projx = zeros(size(coords,1),1);
        projx(1) = coords(1,1);
        for j = 2:size(coords,1);
            D = pdist_euclid([coords(1,:);coords(j,:)]);
            projx(j) = projx(1) + D*cosd(theta);
        end
        projy = [coords(1,2);coords(1,2)+a(i)];
        Gf.(fieldname) = tensorGrid(projx,projy);
        Gf.(fieldname).nodes.coords = new_coords;
        Gf.(fieldname) = computeGeometry(Gf.(fieldname));
    end
    Gf.(fieldname).nodes.start = F(i).nodes.start;
    Gf.(fieldname).cells.start = F(i).cells.start;
    Gf.(fieldname).faces.start = face_start;
    face_start = face_start + Gf.(fieldname).faces.num;
end
G.FracGrid = Gf;
if ~isempty(varargin) && strcmp(varargin{1},'plot')
    figure
    for i = 1:numel(fieldnames(Gf))
        plotGrid(Gf.(['Frac',num2str(i)]),'FaceColor','w');
        hold on
    end
    hold off
end
return


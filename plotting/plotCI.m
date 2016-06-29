function varargout = plotCI(G,fracture)
% plotCI plots conductivity index as matrix cell data for all fracture
% lines 

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


linenum = 1:numel(fracture.lines);
CI = zeros(G.cells.num,numel(fracture.lines));
for i = 1:numel(fracture.lines)
    cells = fracture.lines(linenum(i)).cells;
    for j = 1:numel(cells)
        temp = G.cells.fracture.CI{cells(j),1};
        ls = [G.cells.fracture.line_num{cells(j),1},G.cells.fracture.dfm_line_num{cells(j),1}];
        CI(cells(j),i) = temp(ls==linenum(i));
    end
end
figure
h = plotGrid(G,'FaceColor','none','EdgeAlpha',0.1); hold on
plotToolbar(G,CI,'EdgeAlpha',0.1);
colormap([1 1 1;jet(25)]); colorbar 
ax = gca;
for i = 1:numel(fracture.lines)
    x = [fracture.lines(linenum(i)).endp(1);fracture.lines(linenum(i)).endp(3)];
    y = [fracture.lines(linenum(i)).endp(2);fracture.lines(linenum(i)).endp(4)];
    line(x,y,'Color','k','Parent',ax,'LineWidth',1.5);
end
hold off
title('Fracture-matrix conductivity index per Line');
axis equal tight
if nargout > 0, varargout{1} = h; end
return
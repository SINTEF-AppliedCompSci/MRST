function plotFractureCoarseGrid2D(G,p,F,varargin)
% plotFractureCoarseGrid(G,p,F) is a utility function for plotting, only
% designed to plot the fracture coarse grid with an underlying matrix
% coarse grid. If nargin>3, the script will output global coarse cell
% numbers for the fracture coarse cells beside them. This function uses
% rand() to generate colours for plotting fracture coarse cells.

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


Gm = G.Matrix;
%
frac_cells = transpose(Gm.cells.num+1:G.cells.num);
pfracs = p(frac_cells); % partition number for fracture cells
%
colors = rand(numel(unique(pfracs)),3);
cmap = ones(G.cells.num,3);
if max(p(1:Gm.cells.num))>=max(pfracs)
    for i = 1:numel(frac_cells)
        cmap(frac_cells(i),:) = colors(p(frac_cells(i))==unique(pfracs),:);
    end
else
    cmap(frac_cells,:) = colors(p(frac_cells) - min(pfracs) + 1,:);
end
plotGrid(Gm,'EdgeAlpha',0.04,'FaceColor','none');
outlineCoarseGrid(Gm,p(1:Gm.cells.num),'k');
hold on
for i = 1:numel(F)
    xx = F(i).nodes.coords(:,1);
    yy = F(i).nodes.coords(:,2);
    for j = 1:numel(xx)-1
        clr = cmap(F(i).cells.start+j-1,:);
        plot(xx(j:j+1).',yy(j:j+1).','o-','Color',clr,'LineWidth',2.5,...
            'MarkerEdgeColor',clr,'MarkerFaceColor',clr,'MarkerSize',3);
    end
    if nargin>3 % Plot frac_cell numbers alongside
        xc = 0.5*(xx(1:end-1)+xx(2:end));
        yc = 0.5*(yy(1:end-1)+yy(2:end));
        for j = 1:numel(xc)
            text(xc(j),yc(j),num2str(F(i).cells.start+j-1));
        end
    end
end
dofm = max(p(1:G.Matrix.cells.num));
title({['Fracture coarse grid with ',num2str(numel(unique(pfracs))),' DOF.'],...
    ['Matrix coarse grid with ',num2str(dofm),' DOF.']},...
    'FontSize',15,'FontWeight','bold');
axis equal off tight
hold off
return
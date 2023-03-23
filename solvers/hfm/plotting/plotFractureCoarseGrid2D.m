function varargout = plotFractureCoarseGrid2D(G, p, F, varargin)
% plotFractureCoarseGrid(G, p, F) plots the fracture and matrix coarse
% grids for a 2D domain. This function uses rand() to generate colours for
% plotting fracture coarse cells.
%
% SYNOPSIS:
%       plotFractureCoarseGrid2D(G, p, F)
%       plotFractureCoarseGrid2D(G, p, F, 'pn1', pv1, ...)
%   h = plotFractureCoarseGrid2D(...)
%
% REQUIRED PARAMETERS:
%
%   G  - Grid data structure with fractures as assembled by
%        assembleGlobalGrid.
%
%   p  - Partition vector for the coarse grid.
%
%   F  - see assembleFracNodes2D.
%
% OPTIONAL PARAMETERS:
%
%   showNumbering - If true, plots the fracture coarse block numbers.
%
% RETURNS:
%   h - Handle to polygonal patch structure as defined by function
%       plotFaces.  OPTIONAL.
%
% SEE ALSO:
%   assembleFracNodes2D, outlineCoarseGrid

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

opt = struct('showNumbering', false);

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
h = plotGrid(Gm,'EdgeAlpha',0.04,'FaceColor','none');
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
    if opt.showNumbering % Plot frac_cell numbers alongside
        xc = 0.5*(xx(1:end-1)+xx(2:end));
        yc = 0.5*(yy(1:end-1)+yy(2:end));
        for j = 1:numel(xc)
            text(xc(j),yc(j),num2str(F(i).cells.start+j-1));
        end
    end
end
dofm = max(p(1:G.Matrix.cells.num));
title({['Coarse grid (fracture: ',num2str(numel(unique(pfracs))),' DOFs,'],...
    ['matrix: ',num2str(dofm),' DOFs)']});
axis equal off tight
hold off

if nargout > 0, varargout{1} = h; end

return
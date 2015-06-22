function varargout = plotFractureNodes2D(G,F,fracture,varargin)
% plotFractureNodes2D(G,F,fracture,varargin) plots fracture nodes defined
% in F as returned by assembleFracNodes2D. If varargin{1} = true, it also
% plots global cell numbers beside fractures.

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


dc = distinguishable_colors(numel(fracture.lines),{'w'});
markers = {'o','x','+','*','s','d','v','^','<','>','p','h','.'};
plotGrid(G,'FaceColor','w','EdgeAlpha',0.05);

hold on
count = 1;
for i = 1:numel(F)
    xx = F(i).nodes.coords(:,1);
    yy = F(i).nodes.coords(:,2);
    plot(xx,yy,'Color',dc(i,:),'Marker',markers{count});
    if nargin>3
        xc = 0.5*(xx(1:end-1)+xx(2:end));
        yc = 0.5*(yy(1:end-1)+yy(2:end));
        for j = 1:numel(xc)
            text(xc(j),yc(j),num2str(F(i).cells.start+j-1));
        end
    end
    count = count+1;
    if count == numel(markers)+1
        count = 1;
    end
end
title('Fracture Grid','FontSize',15,'FontWeight','bold');
axis tight
varargout{1} = gcf;
return


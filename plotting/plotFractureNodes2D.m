function varargout = plotFractureNodes2D(G,F,fracture,varargin)
% plotFractureNodes2D plots the 1D fracture grid defined in F as returned
% by assembleFracNodes2D.
%
% SYNOPSIS:
%       plotFractureNodes2D(G,F,fracture)
%       plotFractureNodes2D(G,F,fracture, 'pn1', pv1)
%   h = plotFractureNodes2D(...)
%
% REQUIRED PARAMETERS:
%
%   G        - Matrix grid structure as returned by gridFracture2D.
%
%   F        - Fracture grid structure as returned by gridFracture2D.
%
%   fracture - See gridFracture2D.
%
% OPTIONAL PARAMETERS:
%
%   linewidth     - width of fracture line. Passed as a LineSpec in the
%                   matlab function 'plot'.
%
%   markersize    - Size of markers indicating fracture nodes. Passed as a
%                   LineSpec in the matlab function 'plot'.
%
%   shownumbering - If shownumbering = 1, global cell numbers are also
%                   plotted for the fractures. Might not be readable when
%                   the fracture grid is too fine or there are too many
%                   fractures.
%
% RETURNS:
%   h  - Handle to resulting patch object.  The patch object is added
%        directly to the current AXES object (GCA).
%        OPTIONAL.  Only returned if specifically requested.  If
%        ISEMPTY(cells), then h==-1.
%
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

opt = struct('linewidth'       , 1.5  ,...
             'markersize'      , 3    ,...
             'shownumbering'   , 0    );
opt = merge_options(opt, varargin{:});

dc = rand(numel(fracture.lines),3);
markers = {'o','x','+','*','s','d','v','^','<','>','p','h','.'};
h = plotGrid(G,'FaceColor','w','EdgeAlpha',0.05);

hold on
count = 1;
for i = 1:numel(F)
    xx = F(i).nodes.coords(:,1);
    yy = F(i).nodes.coords(:,2);
    hold on; 
    plot(xx,yy,'Color',dc(i,:),'Marker',markers{count},...
        'MarkerFaceColor', 'k', 'markersize',opt.markersize,...
        'linewidth',opt.linewidth);
    if opt.shownumbering == 1
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
axis tight
if nargout > 0, varargout{1} = h; end
return


function plotFractures(g,cells,data)
%% SYNOPSIS:
%       plotFractures(g)
%       plotFractures(g, data)
%   h = plotFractures(...)
%
% Plot 1D fractures in a 2D grid.
%
% The function plotFaces does not work satisfactory for 2D grids, thus the
% present function can be used instead to visualize fractures.
%
% For 3D grids, plotFaces is applied.
%
% If the grid structure has a field g.faces.fracNodes, these are plotted.
% If not, all faces with g.faces.tags > 0 are ploted instead.
%
% PARAMETERS:
%   g       - Grid data structure.
%   cells   - Index of hybrid cells. If empty, g.cells.hybrid is used
%             instead
%   data    - Data to be plotted. If empty, the fracture geometry is shown
%
% Copyright 2011-2012 University of Bergen
%
% This file is licensed under the GNU General Public License v3.0.

if(nargin)==1 || isempty(cells);
    cells = find(g.cells.hybrid==1);
else
    cells = cells(g.cells.hybrid(cells) == 1);
end

% The linewidth is set to 3, this usually produce reasonable figures. The
% value can be adjusted here
LineWidth = 3;


if g.griddim ==3

    if(nargin)<3
        plotFaces(g,g.cells.tags(cells))
        return
    end
    plotFaces(g,g.cells.tags(cells),data(cells));

else
    if isfield(g.faces,'fracNodes')
        fn = g.faces.fracNodes(g.cells.tags(cells),:);
    else
        fn = getEdgeNodes(g,g.cells.tags(cells));
    end

    if nargin < 3
        patch('Faces',fn,'Vertices',g.nodes.coords,'LineWidth',LineWidth);
        return
    end

    x = [g.nodes.coords(fn(:,1),1)'; g.nodes.coords(fn(:,2),1)'];
    y = [g.nodes.coords(fn(:,1),2)'; g.nodes.coords(fn(:,2),2)'];

    % Assign the face data to both vertices to avoid sorting
    c = repmat(data(cells),1,2)';
    patch(x,y,c,'EdgeColor','flat','FaceColor','none','LineWidth',LineWidth)

end

function edgeNodes = getEdgeNodes(g,faces)
nodePos = g.faces.nodePos(faces);
edgeNodes = [g.faces.nodes(nodePos),g.faces.nodes(nodePos+1)];

function varargout = plotNodeData(G, node_data, varargin)
% Plot data defined on nodes of grid
%
% SYNOPSIS:
%   plotNodeData(G, data)
%
% PARAMETERS:
%   G         - Grid data structure.
%
%   node_data - Data at nodes to be plotted.  
%
%
% KEYWORD ARGUMENTS:
%
%   'Any'    - Additional keyword arguments will be passed directly on to
%              function `patch` meaning all properties supported by `patch`
%              are valid.
%
% RETURNS:
%   h - Handle to resulting PATCH object.  The patch object is added to the
%       current AXES object.
%
% EXAMPLE:
%  G = cartGrid([10, 10]);
%  plotNodeData(G, G.nodes.coords(:, 1));
%
% SEE ALSO:
%   `plotCellData`, `plotGrid`, `newplot`, `patch`, `shading`.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
    opt.cells = [1 : G.cells.num]';
    opt.uu = [];
    [opt, extra] = merge_options(opt, varargin{:});
    if(~isempty(opt.uu))
        G.nodes.coords=G.nodes.coords+opt.uu;
    end
    if (G.griddim == 3)
        [faces, ~] = boundaryFaces(G, opt.cells);
    else
        faces = [1 : G.cells.num]';
    end
    
    % Extract face topology for subset of grid faces and the coordinates for
    % the actual vertices ('verts') present in this topology.
    %
    face_topo  = str2func(sprintf('get_face_topo_%dd', G.griddim));
    [f, verts] = face_topo     (G, faces);
    v          = G.nodes.coords(verts, :);

    if isfield(G.nodes, 'z'), 
        assert (size(v, 2) == 2, ...
                ['Vertex Z coordinate cannot be specified in ', ...
                 'a general 3D grid']);

        v = [v, G.nodes.z(verts)];
    end

    % Massage colour data into a form suitable for 'FaceVertexCData' property.
    %
    % From here, we assume that 'colour' is an m-by-1 (for indexed
    % colouring) or m-by-3 (for 'true colour' colouring) array.
    % Furthermore, 'm' is assumed to be either 1 (meaning all 'faces' should
    % be coloured using a single colour), NUMEL(faces) when the individual
    % faces are coloured separately, or G.nodes.num for interpolated face
    % colouring.
    %
    
    %colour = colour(verts, :);
    fc     = 'interp';
    
    % Build final patch for graphical output (Note: added to GCA).
    
    h = patch('Faces'          , f      , 'Vertices' , v , ...
              'FaceVertexCData', node_data(verts, :) , 'FaceColor', fc, extra{:});
    
    set(get(h, 'Parent'), 'ZDir', 'reverse');
    if nargout > 0
        varargout{1} = h; 
    end
end

%--------------------------------------------------------------------------
% Private helpers follow.
%--------------------------------------------------------------------------

function [f, present] = get_face_topo_2d(G, cells)  %#ok
    eIX = G.cells.facePos;
    nn  = double(diff([G.cells.facePos(cells), ...
                       G.cells.facePos(cells + 1)], [], 2));

    if(isfield(G.cells, 'nodes') )
        cellNodes = G.cells.nodes; 
    else
        cellNodes = getSortedCellNodes(G);
    end

    cn  = double(cellNodes(mcolon(eIX(cells), eIX(cells + 1) - 1), 1));

    m   = numel(cells);
    n   = max(nn);
    f   = nan([n, m]);

    % Extract only those nodes/vertices actually present in the subset of
    % grid faces represented by 'faces'.  Create local numbering for these
    % vertices.
    %
    present           = false([G.nodes.num, 1]);
    present(cn)       = true;

    node_num          = zeros([G.nodes.num, 1]);
    node_num(present) = 1 : sum(double(present));

    off = reshape((0 : m - 1) .* n, [], 1);

    f(mcolon(off + 1, off + nn)) = node_num(cn);

    % PATCH requires that the 'Faces' property be a matrix of size
    % (number of faces)-by-(number of vertices).
    %
    f = f .';
end

%--------------------------------------------------------------------------

function [f, present] = get_face_topo_3d(G, faces)  %#ok
    eIX = G.faces.nodePos;
    nn  = double(diff([G.faces.nodePos(faces), ...
                       G.faces.nodePos(faces + 1)], [], 2));
    fn  = double(G.faces.nodes(mcolon(eIX(faces), eIX(faces + 1) - 1), 1));

    m   = numel(faces);
    n   = max(nn);
    f   = nan([n, m]);

    % Extract only those nodes/vertices actually present in the subset of
    % grid faces represented by 'faces'.  Create local numbering for these
    % vertices.
    %
    present           = false([G.nodes.num, 1]);
    present(fn)       = true;

    node_num          = zeros([G.nodes.num, 1]);
    node_num(present) = 1 : sum(double(present));

    off = reshape((0 : m - 1) .* n, [], 1);

    f(mcolon(off + 1, off + nn)) = node_num(fn);

    tmp = isfinite(f);
    nnode = sum(tmp, 1);
    ind = sub2ind(size(f), nnode, 1 : size(f, 2));
    tmp = repmat(f(ind), size(f, 1), 1);
    f(isnan(f)) = tmp(isnan(f));
    % PATCH requires that the 'Faces' property be a matrix of size
    % (number of faces)-by-(number of vertices).
    %
    f = f .';
end

%--------------------------------------------------------------------------
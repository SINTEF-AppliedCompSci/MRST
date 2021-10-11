function varargout = plotNodeDataDeformed(G, node_data, uu, varargin)
%
% SYNOPSIS:
%       plotNodeDataDeformed(G, node_data, uu, varargin)
%
% DESCRIPTION: Plot node data on a grid which is deformed using a given
% displacement field.
%
% PARAMETERS:
%   G         - Grid data structure.
%   node_data - node data to be plotted  
%
%   'pn'/pv - List of other property name/value pairs.  OPTIONAL.
%             This list will be passed directly on to function PATCH
%             meaning all properties supported by PATCH are valid.
%
% RETURNS:
%   h - Handle to resulting PATCH object.  The patch object is added to the
%       current AXES object.
%
%
% EXAMPLE:
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


    G.nodes.coords = G.nodes.coords + uu;
    if(G.griddim == 3)
        cells = [1 : G.cells.num]';   
        [faces, ~] = boundaryFaces(G, cells);
    else
        faces = [1 : G.cells.num]';
    end
    
    % Extract face topology for subset of grid faces and the coordinates for
    % the actual vertices ('verts') present in this topology.
    %
    face_topo  = str2func(sprintf('get_face_topo_%dd', G.griddim));
    [f, verts] = face_topo(G, faces);
    v          = G.nodes.coords(verts, :);

    if isfield(G.nodes, 'z'), 
        assert (size(v, 2) == 2, ...
                ['Vertex Z coordinate cannot be specified in ', ...
                 'a general 3D grid']);

        v = [v, G.nodes.z(verts)];
    end

    % Massage colour data into a form suitable for 'FaceVertexCData' property.
    %
    % From here, we assume that 'colour' is an m-by-1 (for indexed colouring) or
    % m-by-3 (for 'true colour' colouring) array.  Furthermore, 'm' is assumed to
    % be either 1 (meaning all 'faces' should be coloured using a single colour), 
    % NUMEL(faces) when the individual faces are coloured separately, or
    % G.nodes.num for interpolated face colouring.
    %
    
    % colour = colour(verts, :);
    fc     = 'interp';
    

    if size(f, 1) == size(v, 1)
        % plotFaces uses FaceVertexCData to plot face data. The internal MATLAB
        % patch function decides if the data is given per vertex or face based
        % on a check of lengths. However, if the sizes are exactly *equal* it
        % defaults to vertex data. This can happen when dynamically plotting
        % subsets of fronts. In this case, we add a dummy vertex consisting of
        % NaN which has no influence on the plot, but ensures that the data is
        % interpreted as face data.
        v = [v; NaN(1, size(v, 2))];
    end

    % Build final patch for graphical output (Note : added to GCA).
    
    h = patch('Faces', f, 'Vertices', v, ...
              'FaceVertexCData', node_data(verts, :), ...
              'FaceColor', fc, varargin{:});
    
    set(get(h, 'Parent'), 'ZDir', 'reverse')

    if nargout > 0, varargout{1} = h; end
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

function [plotOutline, varargin] = do_outline_p(varargin)
% Does caller of 'plotFaces' request an outline plot?

    plotOutline = false;
    if numel(varargin) > 0, 
        i = 1 + isnumeric(varargin{1});
        v = find(strcmpi(varargin(i : 2 : end), 'outline'));
        if ~isempty(v), 
            % Find argument following last 'outline'
            plotOutline = varargin{i + 2*v(end) - 1};

            % Remove pairs of 'outline'/value from varargin.
            varargin([i + 2*v-2, i + 2*v - 1]) = [];
        end
    end
end

%--------------------------------------------------------------------------

function rgb = get_rgb(colour)
    switch lower(colour), 
      case {'y', 'yellow' }, rgb = [1, 1, 0];
      case {'m', 'magenta'}, rgb = [1, 0, 1];
      case {'c', 'cyan'   }, rgb = [0, 1, 1];
      case {'r', 'red'    }, rgb = [1, 0, 0];
      case {'g', 'green'  }, rgb = [0, 1, 0];
      case {'b', 'blue'   }, rgb = [0, 0, 1];
      case {'w', 'white'  }, rgb = [1, 1, 1];
      case {'k', 'black'  }, rgb = [0, 0, 0];
      otherwise            , rgb = [0, 0, 1]; % Unknown colour -> 'blue'.
    end
end

%--------------------------------------------------------------------------

function pts = findFaceOutline(g, faces)
    assert (size(g.nodes.coords, 2) == 3);
    % Find points on border of collection of grid faces.
    if numel(faces) == 0, pts = zeros(0, 3); return; end

    cellNo = rldecode(1 : g.cells.num, diff(g.cells.facePos), 2) .';
    sgn    = 2*(cellNo == g.faces.neighbors(g.cells.faces(:, 1), 1)) - 1;
    faceNo = rldecode(1 : g.faces.num, diff(g.faces.nodePos), 2) .';

    fi     = false(g.faces.num, 1); fi(faces) = true;
    fn     = g.faces.nodes(fi(faceNo));

    nodeflag     = false([g.nodes.num, 1]);
    nodeflag(fn) = true;

    fe = faceEdges(g);
    fe(sgn(faceNo)<0, :) = fe(sgn(faceNo)<0, [2, 1]);

    fe = fe (fi(faceNo), :);
    nodes = double(fe(all(nodeflag(fe), 2), :));

    if numel(nodes) > 0, 
        % Remove edges which appear more than once.  These edges are
        % considered internal in the collection of faces.  The remaining
        % edges are on the outer boundary.
        [nodes, n] = rlencode(sortrows(sort(nodes, 2)));
        nodes(n>1, :) = [];
    end

    pts = nan(size(nodes, 1)*3, 3);
    pts (1 : 3 : end, :) = g.nodes.coords(nodes(:, 1), :);
    pts (2 : 3 : end, :) = g.nodes.coords(nodes(:, 2), :);
end

%--------------------------------------------------------------------------

function fe = faceEdges(g)
    fe = [g.faces.nodes, g.faces.nodes([2 : end, 1])];
    fe(g.faces.nodePos(2 : end) - 1, 2) = fe(g.faces.nodePos(1 : end-1), 1);
end

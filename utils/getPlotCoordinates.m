function coords = getPlotCoordinates(G, varargin)
    % Get coordinates for plotting higher-order dG variables

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('n'     , 1000  , ...
                 'plot1d', false);
    
    [opt, ~] = merge_options(opt, varargin{:});

    n0 = sortNodes(G);
    
    if opt.plot1d
        xmax = max(G.nodes.coords, [], 1);
        xmin = min(G.nodes.coords, [], 1);
        x      = linspace(xmin(1), xmax(1), opt.n)';
        y      = repmat((xmin(2) + xmax(2))/2, opt.n, 1);
        points = [x,y];
        cells  = nan(size(points, 1), 1);
        for i = 1:G.cells.num
            xv = G.nodes.coords(G.cells.nodes(G.cells.nodePos(i):G.cells.nodePos(i+1)-1),:);
            ix = inpolygon(points(:,1), points(:,2), xv(:,1), xv(:,2));
            cells(ix) = i;
        end
        faces = [];
    else
    
        faces = find(all(G.faces.neighbors > 0, 2));
        n1 = G.faces.nodes(mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces+1)-1),1);
        n1 = repmat(reshape(n1,2,[]), 2, 1);
        n1 = reshape(n1([1,2,4,3], :), [], 1);
        n = [n0; n1];
        x = G.nodes.coords(n,:);
        cells0 = rldecode((1:G.cells.num)', diff(G.cells.nodePos), 1);


        cells1 = repmat(G.faces.neighbors(faces,:)', 2, 1);
        cells1 = reshape(cells1([1,3,2,4], :), [], 1);
        cells = [cells0; cells1];

        f = (1:size(x,1))';
        ii = [cells0; rldecode((1:numel(faces))' + max(cells0), 4, 1)];
        jj = [mcolon(ones(G.cells.num,1), diff(G.cells.nodePos)), repmat(1:4, 1, numel(faces))]';
        faces = full(sparse(ii, jj, f));
        faces(faces == 0) = nan;
        
    end
    coords = struct('points', x, 'cells', cells, 'faces', faces);

end
    
function n = sortNodes(G)

    f = G.cells.faces(:,1);
    n = G.faces.nodes(mcolon(G.faces.nodePos(f),G.faces.nodePos(f+1)-1));
    s = G.faces.neighbors(f,1) ~= rldecode((1:G.cells.num)', diff(G.cells.facePos),1);

    n = reshape(n, 2, []);
    n(:,s) = n([2,1], s);
    n = n(:);
    n = n(1:2:end);

end

function [cellCenters, cellFaceCenters, cellDims] = ...
      computeCpGeometry(G, grdecl, varargin)
%Undocumented Utility Function

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

opt = struct('on_all_cpc', false);
opt = merge_options(opt, varargin{:});

[fc, cc] = calc_centroids(grdecl);

if ~opt.on_all_cpc
   cellNo      = gridCellNo(G);
   cartCellNo  = G.cells.indexMap(cellNo);
   fType       = G.cells.faces(:, 2);
   cellCenters = cc(G.cells.indexMap, :);
else
   nc = prod(grdecl.cartDims);
   cartCellNo  = rldecode(1:nc, 6, 2) .';
   fType       = repmat((1:6).', [nc, 1]);
   cellCenters = cc;
end

pick_face_centres = @(i) ...
   bsxfun(@times, double(fType == i), fc{i}(cartCellNo, :));

cellFaceCenters  =  ...
   pick_face_centres(1) + pick_face_centres(2) + ...
   pick_face_centres(3) + pick_face_centres(4) + ...
   pick_face_centres(5) + pick_face_centres(6);

if nargout > 2
    vnorm = create_vector_norm_function();

    cellDims = [              ...
       vnorm(fc{2} - fc{1}) , ...
       vnorm(fc{4} - fc{3}) , ...
       vnorm(fc{6} - fc{5}) ];
end
end

%--------------------------------------------------------------------------

function [fc, cc] = calc_centroids(grdecl)
   p = corner_points(grdecl);

   make_face_centre = ...
      @(nodeID) sum(p(:, :, nodeID), 3) ./ numel(nodeID);

   fc{1} = make_face_centre([1, 3, 5, 7]);
   fc{2} = make_face_centre([2, 4, 6, 8]);

   fc{3} = make_face_centre([1, 2, 5, 6]);
   fc{4} = make_face_centre([3, 4, 7, 8]);

   fc{5} = make_face_centre([1, 2, 3, 4]);
   fc{6} = make_face_centre([5, 6, 7, 8]);

   cc = (fc{1} + fc{2}) ./ 2;
end

%--------------------------------------------------------------------------

function p = corner_points(grdecl)
   [x, y, z] = buildCornerPtNodes(grdecl);

   [i, j, k] = ndgrid(1 : 2);
   ijk = [i(:), j(:), k(:)];

   point_coord = @(id) ...
      [coord(ijk(id,:), x), coord(ijk(id,:), y), coord(ijk(id,:), z)];

   p = arrayfun(point_coord, 1 : 8, 'UniformOutput', false);
   p = cat(3, p{:});
end

%--------------------------------------------------------------------------

function c = coord(ijk, v)
    c = reshape(v(ijk(1) : 2 : end, ...
                  ijk(2) : 2 : end, ...
                  ijk(3) : 2 : end), [], 1);
end

%--------------------------------------------------------------------------

function vnorm = create_vector_norm_function()
    persistent has_vecnorm
    if isempty(has_vecnorm)
        has_vecnorm = exist('vecnorm', 'file');
    end
    if has_vecnorm
        vnorm = @(x) vecnorm(x, 2, 2);
    else
        % Fall-back for 'vecnorm'; introduced in MATLAB 9.3.0 (R2017b)
        vnorm = @(x) sqrt(sum(x.^2, 2));
    end
end

function bc = psideh(bc, G, side, fluid, varargin)
%Impose hydrostatic pressure boundary condition on global side.
%
% SYNOPSIS:
%   bc = psideh(bc, G, side, fluid)
%   bc = psideh(bc, G, side, fluid, 'pn', pv)
%   bc = psideh(bc, G, side, fluid, I1, I2)
%   bc = psideh(bc, G, side, fluid, I1, I2, 'pn', pv)
%
% PARAMETERS:
%   bc     - Boundary condition structure as defined by function 'addBC'.
%
%   G      - Grid structure as described by `grid_structure`.  Currently
%            restricted to grids produced by functions `cartGrid` and
%            `tensorGrid` and other grids that add cardinal directions to
%            `G.cells.faces(:, 2) in the same format.
%
%   side   - Global side from which to extract face indices.  String.  Must
%            (case insensitively) match one of six alias groups:
%
%               1) `{'West' , 'XMin', 'Left'  }`
%               2) `{'East' , 'XMax', 'Right' }`
%               3) `{'South', 'YMin', 'Back'  }`
%               4) `{'North', 'YMax', 'Front' }`
%               5) `{'Upper', 'ZMin', 'Top'   }`
%               6) `{'Lower', 'ZMax', 'Bottom'}`
%
%            These groups correspond to the cardinal directions mentioned
%            as the first alternative in each group.
%
%   fluid  - Fluid object.
%
%   I1,I2  - Cell index ranges for local (in-plane) axes one and two,
%            respectively.  An empty index range ([]) is interpreted as
%            covering the entire corresponding local axis of 'side' in the
%            grid 'G'.  The local axes on a 'side' in 'G' are ordered
%            according to 'X' before 'Y', and 'Y' before 'Z'.
%
% OPTIONAL PARAMETERS:
%   sat    - Fluid composition of fluid outside of the reservoir.
%            An m array of fluid phase saturations. If m=3, 'sat'
%            are interpreted as 1 <-> Aqua, 2 <-> Liquid, 3 <-> Vapor.
%
%            This field is to side the density of the outside fluid and to
%            set the saturation of incoming fluid in a transport solver
%
%            Default value: sat = 0 (assume single-phase flow).
%
%   range  - Restricts the search for outer faces to a subset of the cells
%            in the direction perpendicular to that of the face. Example:
%            if side='LEFT', one will only search for outer faces in the
%            cells with logical indexes [:,range,:].
%            Default value: range = [] (do not restrict search)
%
%   ref_depth -
%            Reference depth for pressure. Default is 0.
%
%   ref_pressure -
%            Reference pressure. Default is 0
%
% RETURNS:
%   bc     - Updated boundary condition structure.
%
% EXAMPLES:
%   simpleBC, simpleSRCandBC.
%
% SEE ALSO:
%   `pside`, `fluxside`, `addBC`, `solveIncompFlow`, `grid_structure`.

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

mrstNargInCheck(4, 14, nargin);

if ~isfield(G, 'cartDims')
   error('psideh:NotImplemented', ...
         'psideh is not implemented for this grid type');
end

if nargin==4 || ischar(varargin{1})
   [I1, I2] = deal([]);
else
   I1 = varargin{1};
   I2 = varargin{2};
   varargin = varargin(3:end);
end

opt = struct('sat', 0, 'range', [], ...
             'ref_position', zeros([1, G.griddim]), ...
             'ref_pressure', 0.0);
opt = merge_options(opt, varargin{:});
sat = opt.sat;

ix = boundaryFaceIndices(G, side, I1, I2, opt.range, mfilename());

assert(size(sat,1) == 1, ...
       'Function ''%s'' is only supported with uniform saturations', ...
       mfilename);

[mu, rho] = fluid.properties(struct('s', sat));

kr    = fluid.relperm(sat);
mob   = bsxfun(@rdivide, kr, mu);
omega = dot(mob, rho) ./ sum(mob, 2);

dx = bsxfun(@minus, G.faces.centroids(ix, :), opt.ref_position);
dp = omega .* (dx * reshape(gravity, [], 1));
pressure = opt.ref_pressure + dp;

bc = addBC(bc, ix, 'pressure', pressure, 'sat', sat);

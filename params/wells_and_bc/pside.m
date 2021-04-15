function bc = pside(bc, G, side, pressure, varargin)
%Impose pressure boundary condition on global side.
%
% SYNOPSIS:
%   bc = pside(bc, G, side, p)
%   bc = pside(bc, G, side, p, 'pn', pv)
%   bc = pside(bc, G, side, p, I1, I2)
%   bc = pside(bc, G, side, p, I1, I2, 'pn', pv)
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
%   p      - Pressure value, in units of Pascal, to be applied to the face.
%            Either a scalar or a vector of `numel(I1)*numel(I2)` values.
%
%   I1,I2  - Cell index ranges for local (in-plane) axes one and two,
%            respectively.  An empty index range ([]) is interpreted as
%            covering the entire corresponding local axis of 'side' in the
%            grid 'G'.  The local axes on a 'side' in 'G' are ordered
%            according to 'X' before 'Y', and 'Y' before 'Z'.
%
% OPTIONAL PARAMETERS:
%   sat    - Volumetric composition of fluid injected across inflow faces.
%            An n-by-m array of fluid compositions with 'n' being the
%            number of faces in 'faces' and for m=3, the columns
%            interpreted as 1 <-> Aqua, 2 <-> Liquid, 3 <-> Vapor.
%            The fractions should sum up to one, i.e. have row-sum of
%            unity. If a row vector is specified, this vector is used for
%            all faces in the definition.
%
%            This field is for the benefit of transport solvers such as
%            'implicitTransport' and will be ignored for outflow faces.
%
%            Default value: `sat = 1` (assume single-phase flow).
%
%   range  - Restricts the search for outer faces to a subset of the cells
%            in the direction perpendicular to that of the face. Example:
%            if side='LEFT', one will only search for outer faces in the
%            cells with logical indexes [range,:,:].
%            Default value: `range = []` (do not restrict search).
%
% RETURNS:
%   bc     - Updated boundary condition structure.
%
% EXAMPLE:
%   See simpleBC, simpleSRCandBC.
%
% SEE ALSO:
%   `fluxside`, `addBC`, `grid_structure`.

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


if ~isfield(G, 'cartDims')
   error(msgid('NotImplemented'), ...
         'PSIDE is not implemented for this grid type.');
end

mrstNargInCheck(4, 10, nargin);

if nargin == 4 || ischar(varargin{1})
   % pside(bc, G, side, p, ['pn1', pv1, ...]).  Use entire face.
   I1 = []; I2 = [];
else
   % pside(bc, G, side, p, I1, I2, ['pn1', pv1, ...]).
   I1 = varargin{1}; I2 = varargin{2};
   varargin = varargin(3 : end);
end

opt = struct('sat', 1, 'range', []);
opt = merge_options(opt, varargin{:});
sat = opt.sat;

ix = boundaryFaceIndices(G, side, I1, I2, opt.range, mfilename());

assert (any(numel(pressure) == [1, numel(ix)]));
assert (isempty(sat) || any(size(sat,1) == [1, numel(ix)]));

if size(sat,1)     == 1, sat      = sat(ones([numel(ix), 1]), :); end
if numel(pressure) == 1, pressure = pressure(ones([numel(ix), 1])); end

bc = addBC(bc, ix, 'pressure', pressure, 'sat', sat);

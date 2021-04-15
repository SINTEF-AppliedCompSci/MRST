function bc = addBC(bc, f, t, v, varargin)
%Add boundary condition to (new or existing) BC object
%
% There can only be a single boundary condition per face in the grid.
% Solvers assume boundary conditions are given on the boundary; conditions
% in the interior of the domain yield unpredictable results and is not
% officially supported. Faces with no boundary conditions are generally
% interpreted as no-flux boundary conditions for all phases.
%
% SYNOPSIS:
%   bc = addBC(bc, faces, type, values)
%   bc = addBC(bc, faces, type, values, 'pn1', pv1, ...)
%
% PARAMETERS:
%   bc     - Boundary condition structure from a prior call to `addBC`
%            which will be updated on output or an empty array (`bc==[]`)
%            in which case a new boundary condition structure is created.
%
%   faces  - Global faces in external model for which this boundary
%            condition should be applied.
%
%   type   - Type of boundary condition.  Supported values are 'pressure'
%            and 'flux', or cell array of such strings.
%
%   values - Boundary condition value.  Interpreted as a pressure value (in
%            units of 'Pa') when `type=='pressure'` and as a flux value (in
%            units of 'm^3/s') when `type=='flux'`.  One scalar value for
%            each face in 'faces'. If a single value is given, it will be
%            repeated for all faces.
%            
%
%            Note: If `type=='flux'`, the values in 'values' are interpreted
%            as injection fluxes (into the reservoir).  Specifically, a
%            positive value is interpreted as an injection flux.  To
%            specify an extraction flux (i.e., flux out of the reservoir),
%            the caller should provide a negative value in 'values'.
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
% NOTE:
%   For convenience, values and sat may contain a single value.  This value
%   is then used for all faces specified in the call.
%
% RETURNS:
%   bc - New or updated boundary condition structure having the following
%        fields:
%          - face:  External faces for which explicit BCs are provided.
%          - type:  Cell array of strings denoting type of BC.
%          - value: Boundary condition values for all faces in 'face'.
%          - sat: Fluid composition of fluids passing through inflow faces.
%
% SEE ALSO:
%   `pside`, `fluxside`, `addSource`, `addWell`, `solveIncompFlow`, `grid_structure`.

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


if isempty(f)
   warning('MRST:addBC', 'Empty list of boundary faces.');
   return;
end

opt = struct('sat', 1);
opt = merge_options(opt, varargin{:});
s   = opt.sat;

if isempty(bc)
   bc = struct('face', [], 'type', {{}}, 'value', [], 'sat', []);
end

% Validate boundary condition type.
% Convert to (one-element) cell array if valid.
%
if ischar(t)
   t = repmat({ t }, [1, numel(f)]);
else
   assert(iscellstr(t) && numel(t) == numel(f), ...
      ['Boundary condition type should be either a string or a cell', ...
       'array of strings']);
    t = reshape(t, 1, []);
end

t = lower(t);

assert (all(strcmpi(t, 'pressure') | strcmpi(t, 'flux')), ...
        'Boundary condition type should be either ''pressure'' of ''flux''');

% Validate saturation input.
%  - If ISEMPTY(bc.sat), then 's' may be any numeric array (including []).
%  - Otherwise, 's' must be numeric and have the same number of columns as
%    the existing 'bc.sat' array.
%
assert (isnumeric(s));
assert (isempty(bc.sat) || (size(bc.sat,2) == size(s,2)),...
    'Number of columns in sat field does not match pre-existing sat field in bc.');

% Verify that boundary condition is not already set
bc_given = false(max([f(:); bc.face]), 1);
bc_given(bc.face) = true;
assert(~any(bc_given(f)), ...
   'New boundary condition overlaps with conditions already in bc-struct.');

nf = numel(f);
% Expand single-element saturations and BC values to cover all faces.
%
if size(s,1) == 1
    s = repmat(s, nf, 1);
end
if numel(v)  == 1
    v = repmat(v, nf, 1);
end

% Verify that v and s are same length as faces (or s empty).
assert (numel(v)      ==     nf  );
assert (any(size(s,1) == [0, nf]));

% Boundary conditions structurally verified.  Append to existing structure.
bc.face  = [bc.face ; f(:)];
bc.type  = [bc.type , t];
bc.value = [bc.value; v(:)];
bc.sat   = [bc.sat  ; s];

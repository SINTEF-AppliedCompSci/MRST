function src = addSource(src, c, r, varargin)
%Add an explicit source to (new or existing) source object.
%
% There can only be a single net source term per cell in the grid.  This
% is now enforced.
%
% SYNOPSIS:
%   src = addSource(src, cells, values)
%   src = addSource(src, cells, values, 'pn1', pv1)
%
% PARAMETERS:
%   src    - Source structure from a prior call to 'addSource' which will
%            be updated on output or an empty array (`src==[]`) in which 
%            case a new source structure is created.
%
%   cells  - Indices in external model for which this source should be
%            applied.
%
%   values - Strength of source.  One scalar value for each cell in
%            'cells'.  Note that the values in 'values' are interpreted as
%            flux rates (typically in units of m^3/Day) rather than as flux
%            density rates (which must be integrated over the cell volumes
%            in order to obtain flux rates).  Specifically, the mimetic
%            pressure solvers do not integrate these values.
%
%            In the special case that a single value is provided, it will
%            be assumed valid for all cells in the input.
%
% OPTIONAL PARAMETERS:
%   sat    - Fluid composition of injected fluid in injection cells.
%            An n-by-m array of fluid compositions with 'n' being the
%            number of cells in 'cells' and 'm' the number of components in
%            the saturation. For m=3, the columns interpreted as
%              1 <-> Aqua, 2 <-> Liquid, 3 <-> Vapor.
%
%            This field is for the benefit of transport solvers such as
%            `implicitTransport` and will be ignored for production source
%            cells (i.e. when values < 0).
%
%            Default value: sat = [] (assume single-phase flow).
%
%            As a special case, if `size(sat,1)==1`, then the saturation
%            value will be repeated for all affected cells defined by the
%            'cells' parameter.
%
% RETURNS:
%   src - New or updated source structure having the following fields:
%          - cell: Cells for which explicit sources are provided
%          - rate: Rates or values of these explicit sources
%          - sat:  Fluid composition of injected fluids in injection cells.
%              
%
% EXAMPLE:
%   simpleSRCandBC.m
%
% SEE ALSO:
%   `addWell`, `addBC`, `solveIncompFlow`.

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


opt = struct('sat', []);
opt = merge_options(opt, varargin{:});
s   = opt.sat;

nc = numel(c);
if ~isempty(s) && size(s,1) == 1
   % Support single sat in all affected cells.
   s = repmat(s, nc, 1);
end

if isempty(src)
   src = struct('cell', [], 'rate', [], 'sat', []);
end

if numel(r) == 1
    r = repmat(r, nc, 1);
end

assert (numel(c) == numel(r), ['The number of rates should equal the number', ...
                    ' of cells on input, or be a single value for all cells.']);
assert ((size(s,2) == size(src.sat,2)) || (size(src.sat,2) == 0));

% Verify that cell source term is not already set
src_given = false(max([c(:); src.cell]), 1);
src_given(src.cell) = true;
assert(~any(src_given(c(:))), ...
   'New source terms in conflict with source terms stored in src-struct.');


src.cell = [src.cell; c(:)];
src.rate = [src.rate; r(:)];
src.sat  = [src.sat ; s];

function ix = boundaryFaceIndices(G, direction, i1, i2, i3, caller)
%Retrieve face indices belonging to subset of global outer faces.
%
% SYNOPSIS:
%   ix = boundaryFaceIndices(G, side, i1, i2, i3)
%
% PARAMETERS:
%   G     - Grid data structure.
%
%   side  - Global side from which to extract face indices.  String.  Must
%           (case insensitively) match one of six alias groups:
%
%              1) {'West' , 'XMin', 'Left'  }
%              2) {'East' , 'XMax', 'Right' }
%              3) {'South', 'YMin', 'Back'  }
%              4) {'North', 'YMax', 'Front' }
%              5) {'Upper', 'ZMin', 'Top'   }
%              6) {'Lower', 'ZMax', 'Bottom'}
%
%           These groups correspond to the cardinal directions mentioned as
%           the first alternative in each group.
%
% OPTIONAL PARAMETERS:
%
%   i1,i2 - Index ranges for local (in-plane) axes one and two,
%           respectively.  An empty index range ([]) is interpreted as
%           covering the entire corresponding local axis of 'side' in the
%           grid 'G'.  The local axes on a 'side' in 'G' are ordered
%           according to 'X' before 'Y', and 'Y' before 'Z'.
%
%   i3    - Index range for global axis perpendicular to 'side'.  The
%           primary purpose of this parameter is to exclude faces *within*
%           a reservoir from being added to the return value 'ix'.  Such
%           faces typically occur in faulted reservoirs where a given face
%           may be considered external by virtue of being connected to a
%           single reservoir cell only.
%
% RETURNS:
%   ix    - Required face indices.
%
%
% SEE ALSO:
%   `fluxside`, `pside`.

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

if nargin < 3
    i1 = [];
else
    i1 = reshape(i1, [], 1);
end
if nargin < 4
    i2 = [];
else
    i2 = reshape(i2, [], 1);
end
if nargin < 5
    i3 = [];
else
    i3 = reshape(i3, [], 1);
end

if nargin < 6
    % Undocumented input - used to preserve the ID in errors/warnings from
    % pside/fluxside.
    caller = mfilename();
end

%--------------------------------------------------------------------------
% Extract all faces of cells within the given subset.
%
[cells, ft, isOutF] = bdryCellsSubset(G, direction, i1, i2, i3, caller);

fIX   = G.cells.facePos;
hfIX  = mcolon(fIX(cells), fIX(cells + 1)-1);
faces = G.cells.faces(hfIX, 1);
tags  = G.cells.faces(hfIX, 2);

%--------------------------------------------------------------------------
% Extract only those faces which have the required tag.
%
ix = faces(isOutF(faces) & tags == ft);


%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------


function [cells, faceTag, isOutF] = ...
      bdryCellsSubset(G, direction, i1, i2, i3, caller)
% Determine which indices and face tags to look for
switch lower(direction)
   case {'left'  , 'west' , 'xmin'}
      % I == min(I)
      [d1, d2, d3, faceTag] = deal(2, 3, 1, 1);

   case {'right' , 'east' , 'xmax'}
      % I == max(I)
      [d1, d2, d3, faceTag] = deal(2, 3, 1, 2);

   case {'back'  , 'south', 'ymin'}
      % J == min(J)
      [d1, d2, d3, faceTag] = deal(1, 3, 2, 3);

   case {'front' , 'north', 'ymax'}
      % J == max(J)
      [d1, d2, d3, faceTag] = deal(1, 3, 2, 4);

   case {'top'   , 'upper', 'zmin'}
      % K == min(K)
      [d1, d2, d3, faceTag] = deal(1, 2, 3, 5);

   case {'bottom', 'lower', 'zmax'}
      % K == max(K)
      [d1, d2, d3, faceTag] = deal(1, 2, 3, 6);

   otherwise
      error([caller, ':Side:Unknown'],                  ...
            'Boundary side ''%s'' not supported in %s', ...
            direction, caller);
end

% Determine uniqe outer cells (i.e., cells on boundary)
isOutF = any(double(G.faces.neighbors) == 0, 2);
cells  = find(accumarray(sum(G.faces.neighbors(isOutF,:), 2), 1) > 0);

% Determine logical indices of these cells
% Assume we will only be called for logically Cartesian grids for which the
% fields 'G.cartDims' and 'G.cells.indexMap' are present.
%
dims = reshape(G.cartDims, 1, []);
if G.griddim == 1
   dims = [dims(1), 1, 1];

   if faceTag > 2
      error([caller, ':Side:Unsupported'], ...
            ['Boundary side ''%s'' is not defined for one-dimensional ', ...
             'grids in ''%s''.'], direction, caller);
   end
end

if G.griddim == 2
   dims = [dims(1:2), 1];

   if faceTag > 4
      error([caller, ':Side:Unsupported'], ...
            ['Boundary side ''%s'' is not defined for two-dimensional ', ...
             'grids in ''%s''.'], direction, caller);
   end

   if numel(i2) > 1
      error([caller, ':I2:ERANGE'], ...
            ['Two-dimensional boundary faces are incompatible with a\n', ...
             'two-dimensional grid model in ''%s''.\n\n'               , ...
             'Specifically, ''I2'' must contain only a single number.'], ...
            caller);
   end

   if numel(i3) > 0 && any(i3 > 1)
      error([caller, ':I3:ERANGE'], ...
            ['A non-zero cell depth is incompatible with a\n', ...
             'two-dimensional grid model in ''%s''.\n\n'     , ...
             'Specifically, ''I3'' must be empty.'], caller);
   end
end
[cIJK{1:3}] = ind2sub(dims, double(G.cells.indexMap(cells)));

if isempty(i1)
    i1 = 1:dims(d1);
end

if isempty(i2)
    i2 = 1:dims(d2);
end

% Determine whether or not a given cell is within the required subset
%
if any(i1 < 1 | dims(d1) < i1)
   error([caller, ':I1:ERANGE'], ...
          'Cell range ''I1'' outside model in ''%s''.', caller);
end
if any(i2 < 1 | dims(d2) < i2)
   error([caller, ':I2:ERANGE'], ...
          'Cell range ''I2'' outside model in ''%s''.', caller);
end
if numel(i3) > 0 && any(i3 < 1 | dims(d3) < i3)
   error([caller, ':I3:ERANGE'], ...
          'Cell range ''I3'' outside model in ''%s''.', caller);
end

[I{1:2}] = deal(false([G.cells.num, 1]));
I{1}(i1) = true;
I{2}(i2) = true;
inSubSet = I{1}(cIJK{d1}) & I{2}(cIJK{d2});

if ~isempty(i3)
   I{3}     = false([G.cells.num, 1]);   I{3}(i3) = true;
   inSubSet = inSubSet & I{3}(cIJK{d3});
end

% Extract required cell subset
cells = cells(inSubSet);

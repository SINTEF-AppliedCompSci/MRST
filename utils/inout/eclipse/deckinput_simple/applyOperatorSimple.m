function sect = applyOperatorSimple(sect, fid, cartDims, box, kw)
%Apply ECLIPSE/FrontSim operator to input array.
%
% SYNOPSIS:
%   sect = applyOperatorSimple(sect, fid, cartDims, box, kw)
%
% PARAMETERS:
%   sect - Structure array representing the current deck section.
%
%   fid  - Valid file identifier (as obtained using FOPEN) to file
%          containing array operator description.  FTELL(fid) is assumed to
%          be at the start of the description (i.e., after the keyword
%          which prompted the reading of this table).
%
%   cartDims -
%          Cartesian grid dimensions.  A three element vector of integers
%          specifying the number of cells in the 'X', 'Y', and 'Z'
%          directions, respectively.  Typically corresponds to the 'DIMENS'
%          keyword of the input deck.
%
%   box  - Current input box.  Presently assumed to equal 'cartDims'.
%
%   kw   - Operator keyword.  Must be one of
%             ADD      -- Add constant value to existing array.
%             COPY     -- Copy values from existing array to other array.
%             EQUALS   -- Assign constant value to array.
%             MAXVALUE -- Array elements must not exceed specific value.
%             MINVALUE -- Array elements must not be less than given value.
%             MULTIPLY -- Multiply existing array by constant value.
%
% RETURNS:
%   sect - Updated deck section structure array.
%
% SEE ALSO:
%   `readGRID`, `readEDIT`, `readPROPS`, `readREGIONS`, `readSOLUTION`.

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


   assert (any(strcmp(kw, {'ADD', 'COPY', 'EQUALS', ...
                           'MAXVALUE', 'MINVALUE', 'MULTIPLY'})));

   tmpl(1 : 8)       = { '1' };
   tmpl(4 : 2 : end) = arrayfun(@num2str, box, 'UniformOutput', false);

   op = str2func(kw);

   rec = readDefaultedRecord(fid, tmpl);
   while ~isequal(rec, tmpl),
      [idx, m] = box_indices(to_double(rec(3 : end)), cartDims);

      sect = op(sect, rec{1 : 2}, idx, m);

      rec = readDefaultedRecord(fid, tmpl);
   end
end

%--------------------------------------------------------------------------

function sect = ADD(sect, x, val, ix, m)
   % x <- x + val
   assert (isfield(sect, x) && numel(sect.(x)) >= m);
   assert (~any(isnan(sect.(x)(ix))), ...
           'Cannot add value to undefined array ''%s''.', x);

   sect.(x)(ix) = sect.(x)(ix) + to_double(val);
end

%--------------------------------------------------------------------------

function sect = COPY(sect, x, y, ix, m)
   % y <- x
   assert (isfield(sect, x) && numel(sect.(x)) >= m);

   if ~isfield(sect, y),
      % Allocate if not present.
      sect.(y) = zeros([max(numel(sect.(x)), m), 1]);
   end

   sect.(y)(ix) = sect.(x)(ix);
end

%--------------------------------------------------------------------------

function sect = EQUALS(sect, x, val, ix, m)
   % x <- val
   if ~isfield(sect, x),
      % Allocate if not present.
      sect.(x) = nan([m, 1]);
   elseif numel(sect.(x)) < m,
      sect.(x)(numel(sect.(x)) + 1 : m) = nan;
   end

   sect.(x)(ix) = to_double(val);
end

%--------------------------------------------------------------------------

function sect = MAXVALUE(sect, x, val, ix, varargin)
   % x <- min(x, val), i.e., x must not exceed 'val'
   if isfield(sect, x),
      sect.(x)(ix) = min(sect.(x)(ix), to_double(val));
      sect.(x)(ix(isnan(sect.(x)(ix)))) = to_double(val);
   end
end

%--------------------------------------------------------------------------

function sect = MINVALUE(sect, x, val, ix, varargin)
   % x <- max(x, val), i.e., x must not be less than 'val'
   if isfield(sect, x),
      sect.(x)(ix) = max(sect.(x)(ix), to_double(val));
      sect.(x)(ix(isnan(sect.(x)(ix)))) = to_double(val);
   end
end

%--------------------------------------------------------------------------

function sect = MULTIPLY(sect, x, val, ix, m)
   % x <- x * val
   assert (isfield(sect, x) && numel(sect.(x)) >= m);
   assert (~any(isnan(sect.(x)(ix))), ...
           'Cannot multiply undefined array ''%s''.', x);

   sect.(x)(ix) = sect.(x)(ix) .* to_double(val);
end

%--------------------------------------------------------------------------

function [ix, m] = box_indices(i, dims)
   [j{1 : 3}] = ndgrid(i(1) : i(2), i(3) : i(4), i(5) : i(6));

   ix = sub2ind(reshape(dims, 1, []), ...
                reshape(j{1}, [], 1), ...
                reshape(j{2}, [], 1), ...
                reshape(j{3}, [], 1));
   m  = ix(end);  % == MAX(ix) in this case.
end

%--------------------------------------------------------------------------

function v = to_double(v)
   assert (ischar(v) || iscell(v));
   convert = @(s) sscanf(s, '%f');

   if ischar(v),
      v = convert(v);
   else
      v = reshape(cellfun(convert, v), [], 1);
   end
end

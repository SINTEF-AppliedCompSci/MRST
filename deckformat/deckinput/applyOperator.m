function sect = applyOperator(sect, fid, kw)
%Apply ECLIPSE/FrontSim operator to input array.
%
% SYNOPSIS:
%   sect = applyOperator(sect, fid, kw)
%
% PARAMETERS:
%   sect - Structure array representing the current deck section.
%
%   fid  - Valid file identifier (as obtained using FOPEN) to file
%          containing array operator description.  FTELL(fid) is assumed to
%          be at the start of the description (i.e., after the keyword
%          which prompted the reading of this table).
%
%   kw   - Operator keyword.  Must be one of
%             ADD      -- Add constant value to existing array.
%             COPY     -- Copy values from existing array to other array.
%             EQUALS   -- Assign constant value to array.
%             MAXVALUE -- Array elements must not exceed specific value.
%             MINVALUE -- Array elements must not be less than given value.
%             MULTIPLY -- Multiply existing array by constant value.
%
% NOTE:
%   Box limits for the operator 'kw' are managed independently from the
%   limits that affect the property keywords (e.g., PERMX).  The limits are
%   initialised to that of the most recent input box defined by the BOX or
%   ENDBOX keywords--or the entire model if no such box limits have been
%   previously selected in the input deck.
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
                           'MAXVALUE', 'MINVALUE', 'MULTIPLY'})), ...
          ['Keyword ''%s'' not amongst supported operators.\n', ...
           'Must be one of\n\t', ...
           '* ADD\n\t'         , ...
           '* COPY\n\t'        , ...
           '* EQUALS\n\t'      , ...
           '* MAXVALUE\n\t'    , ...
           '* MINVALUE\n\t'    , ...
           '* MULTIPLY.'], kw);

   tmpl(1 : 8) = { '1' };

   dims = defaultBox;
   dims = dims(2 : 2 : end);
   m    = prod(dims);

   op = str2func(kw);

   [rec, tmpl, b] = get_record(fid, tmpl, gridBox);

   while ~isequal(rec, tmpl)
      target = rec{1};
      if ~isfield(sect, target)
         switch kw
            case 'ADD'
                sect.(target) = zeros([m, 1]);
            case 'MULTIPLY'
                sect.(target) = ones([m, 1]);
            case 'EQUALS'
                if regexpi(target, 'MULT\w+')
                    % Multipliers should be initialized as one
                    v = ones([m, 1]);
                else
                    v = nan([m, 1]);
                end
                sect.(target) = v;
            otherwise
               warning(msgid('Default:Unavailable'), ...
                     ['No reasonable initial value for operator ', ...
                      '''%s'' applied to undefined array ''%s''.'], ...
                      kw, target);
         end
      end

      idx = box_indices(dims, b);

      sect = op(sect, rec{1 : 2}, idx, m);

      [rec, tmpl, b] = get_record(fid, tmpl, b);
   end
end

%--------------------------------------------------------------------------

function [rec, tmpl, b] = get_record(fid, tmpl, b)
   tmpl(3 : end) = arrayfun(@int2str, b, 'UniformOutput', false);

   rec = readDefaultedRecord(fid, tmpl);

   b = cellfun(@str2double, rec(3 : end));
end

%--------------------------------------------------------------------------

function ix = box_indices(dims, b)
   [i{1 : 3}] = ndgrid(b(1) : b(2), b(3) : b(4), b(5) : b(6));

   ix = sub2ind(reshape(dims, 1, []), ...
                reshape(i{1}, [], 1), ...
                reshape(i{2}, [], 1), ...
                reshape(i{3}, [], 1));
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
   assert (isfield(sect, x) && numel(sect.(x) >= m));

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

function v = to_double(v)
   assert (ischar(v) || iscell(v));
   convert = @(s) sscanf(s, '%f');

   if ischar(v),
      v = convert(v);
   else
      v = reshape(cellfun(convert, v), [], 1);
   end
end

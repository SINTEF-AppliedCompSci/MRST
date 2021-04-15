function pvto = completePVTO(pvto)
%Fill in missing u-sat data from scaled copies of the FVF and viscosity.
%
% SYNOPSIS:
%   pvto = completePVTO(pvto)
%
% PARAMETERS:
%   pvto - ECLIPSE ``live oil'' PVT table as defined by function
%          'readMisciblePVTTable', possibly lacking undersaturated data for
%          some records.
%
% RETURNS:
%   pvto - Updated ``live oil'' table where all records provide u-sat data.
%          Missing data is provided by scaled copies of the formation
%          volume factor and viscosity of later records.  The scaling is
%          done in such a way that the compressibility and viscosibility of
%          the data is preserved.
%
% NOTE:
%   This is an internal function that is not (and should never be) exposed
%   to user code.
%
% SEE ALSO:
%   `readMisciblePVTTable`, `private/completePVTG`.

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


   for i = 1 : numel(pvto),
      d = diff(pvto{i}.pos);

      if any(d == 1),
         pvto{i} = fill_usat_data(pvto{i}, d);
      end
   end
end

%--------------------------------------------------------------------------
% Private helpers follow.
%--------------------------------------------------------------------------

function t = fill_usat_data(t, d)
   [pos, data, add, ix, n, i] = expand_existing(t, d);
   data(i, :) = scaled_copy(t, d, add, ix, n);

   t.pos  = pos;
   t.data = data;
end

%--------------------------------------------------------------------------

function [pos, data, add, ix, n, i2] = expand_existing(t, d)
   add = d == 1;        % Records missing usat data
   cpy = find(~add);    % Records providing usat data (viable for copying)
   add = find( add);

   % Compute index into 'cpy' (i.e., usat data record number) from which to
   % copy the next record.
   %
   % Algorithm:
   %
   %   1) BSXFUN(@GT, add, cpy) -> NUMEL(add)-by-NUMEL(cpy) array of 0/1
   %   2) SUM(..., 2)           -> Number of 'cpy' entries less than 'add'
   %   3) 1 + ...               -> Index of first viable 'cpy' record
   %   4) RLENCODE(...)         -> Compress/count unique copy rows.
   %
   % 'ix' is index, 'n' is repeat count.
   %
   [ix, n] = rlencode(1 + sum(bsxfun(@gt, add, reshape(cpy, 1, [])), 2));

   % Compute row index for (linearised) copied data (i).  Second output is
   % not used, but needed for 'blockDiagIndex' semantics.  This definition
   % supports multi-row usat copying.
   %
   [i, j] = blockDiagIndex(ix, n);  %#ok

   % Count number of extra rows in each record (one less than number of
   % rows of copied record).
   extra      = zeros([numel(add) + numel(cpy), 1]);
   extra(add) = d(cpy(rldecode(ix, n))) - 1;

   % Compute indirection array (pos2) into expanded data array.  Allocate
   % the expanded array (data).
   pos  = cumsum([1; d + extra]);
   data = zeros([pos(end)-1, size(t.data, 2)]);
   ix   = cpy(ix);

   % Define 'data' row indices into which to fill original data (i1), and
   % copied, usat data (i2).  The 'i1' index is simply the original indices
   % suitably translated for insertion of copied data.  Similarly, d(add)
   % is the offset at which to start filling the copied data.
   %
   i1 = mcolon(pos(1 : end - 1) , pos(1:end-1) + d - 1);
   i2 = mcolon(pos(add) + d(add), pos(add + 1)     - 1);

   % Fill original data into expanded table.
   data(i1, :) = t.data;
end

%--------------------------------------------------------------------------

function data = scaled_copy(t, d, add, ix, n)
   nloc = diff([t.pos(ix), t.pos(ix + 1)], [], 2);

   offset = cumsum([1; nloc(1 : end - 1)]);
   extent = nloc - 1;

   X = t.data(mcolon(t.pos(ix), t.pos(ix + 1) - 1), :);

   X(:, 2) = 1 ./ X(:, 2);      % Interpolate 1/B

   dX = diff(X);

   % Compute compressibilities and viscosibilities (C).
   % Derivatives approximated at left end point.
   %
   I = mcolon(offset, offset + extent - 1);
   C = dX(I, 2:end) ./ X(I, 2:end);
   C = bsxfun(@rdivide, C, dX(I, 1));

   % Preserve C when copying u-sat data.
   %
   [I, J] = blockDiagIndex(extent, n);   %#ok
end


%{
%--------------------------------------------------------------------------

function t = copy_undersat_data(t)
   n0  = diff(t.pos);   % Number of rows in each (orig) record
   add = n0 == 1;       % Records missing usat data
   cpy = find(~add);    % Records providing usat data (viable for copying)
   add = find( add);

   % Compute index into 'cpy' (i.e., usat data record number) from which to
   % copy the next record.
   %
   % Algorithm:
   %
   %   1) BSXFUN(@GT, add, cpy) -> NUMEL(add)-by-NUMEL(cpy) array of 0/1
   %   2) SUM(..., 2)           -> Number of 'cpy' entries less than 'add'
   %   3) 1 + ...               -> Index of first viable 'cpy' record
   %   4) RLENCODE(...)         -> Compress/count unique copy rows.
   %
   % 'ix' is index, 'n' is repeat count.
   %
   [ix, n] = rlencode(1 + sum(bsxfun(@gt, add, reshape(cpy, 1, [])), 2));

   % Compute row index for (linearised) copied data (i).  Second output is
   % not used, but needed for 'blockDiagIndex' semantics.  This definition
   % supports multi-row usat copying.
   %
   [i, j] = blockDiagIndex(ix, n);  %#ok

   % Count number of extra rows in each record (one less than number of
   % rows of copied record).
   extra      = zeros([numel(add) + numel(cpy), 1]);
   extra(add) = n0(cpy(rldecode(ix, n))) - 1;

   % Compute indirection array (pos2) into expanded data array.  Allocate
   % the expanded array (data).
   pos2 = cumsum([1; n0 + extra]);
   data = zeros([pos2(end)-1, size(t.data, 2)]);

   % Extract usat data from the requested records (linearised).
   dcpy = t.data(mcolon(t.pos(cpy(ix)    ) + 1, ...
                        t.pos(cpy(ix) + 1) - 1), : );

   % Define 'data' row indices into which to fill original data (i1), and
   % copied, usat data (i2).  The 'i1' index is simply the original indices
   % suitably translated for insertion of copied data.  Similarly, n0(add)
   % is the offset at which to start filling the copied data.
   %
   i1 = mcolon(pos2(1 : end - 1)  , pos2(1:end-1) + n0 - 1);
   i2 = mcolon(pos2(add) + n0(add), pos2(add + 1)      - 1);

   % Fill original and copied data.
   data(i1, :) = t.data;
   data(i2, :) = dcpy(i,:);

   % Define final table data.
   t.data = data;
   t.pos  = pos2;
end
%}

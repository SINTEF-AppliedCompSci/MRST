function T = readMisciblePVTTable(fid, ntab, ncol, kw)
%Input miscible ECLIPSE/FrontSim PVT table (e.g., PVTO, PVTG).
%
% SYNOPSIS:
%   Table = readMisciblePVTTable(fid, ntab, ncol)
%
% PARAMETERS:
%   fid  - Valid file identifier (as obtained from FOPEN) to file
%          containing the PVT table.  FTELL(fid) is assumed to be at the
%          start of the PVT table (i.e., after the keyword, if any, which
%          prompted the reading of this table).
%
%   ntab - Number of PVT tables to input.  Typically corresponds to the
%          ECLIPSE/FrontSim parameter 'NTPVT' (item two of 'TABDIMS').
%
%   ncol - Number of columns in each of the 'ntab' tables.
%
% RETURNS:
%   Table - An ntab-by-1 cell array of miscible PVT tables.  A miscible PVT
%           table is represented as a data structure (STRUCT) with the
%           following fields:
%
%             - key  -- Primary lookup key for data table interpolation.
%                       This is typically pressure or Rs values.  The key
%                       may be used for binning input quantities before
%                       interpolating the 'data' portion related to each
%                       key.  For binning purposes, 'key' should normally
%                       be sorted in the input source represented by 'fid'.
%
%             - pos  -- Indirection map into 'data'.  Specifically, the
%                       data related to 'key(i)' is found in the submatrix
%
%                          data(pos(i) : pos(i + 1) - 1, :)
%
%             - data -- Table data.  An m-by-ncol array of type DOUBLE.
%                       It is the caller's responsibility to correctly
%                       interpret the individual columns of a specific
%                       table.
%
% NOTE:
%   Miscible PVT tables in ECLIPSE or FrontSim contain an arbitrary number
%   of records, each record terminated by a slash character ('/').  At
%   least one record (the last), and optionally all other records, must
%   specify additional data (typically Rv or pressure values for
%   undersaturated fluid solutions).
%
%   Function 'readMisciblePVTTable' does not verify or enforce sorted keys.
%
%   The return value, 'Table', may be converted to a structure array using
%   the statement:
%
%         Table = [ Table{:} ] .'  % or:  Table = vertcat(Table{:})
%
% SEE ALSO:
%   `readRelPermTable`, `readImmisciblePVTTable`.

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


   T = cell([ntab, 1]);

   for i = 1 : ntab
      t    = create_table(ncol);
      done = false;

      while ~done
         data      = read_record(fid);
         [t, done] = append_data(t, data, ncol);
      end

      T{i} = fill_usat_invlinear(finish_table(t));
   end
end

%--------------------------------------------------------------------------

function data = read_record(fid)
   data = readRecordString(fid);

   if isempty(data) || all(isspace(data))
      data = [];
   else
      data = sscanf(data, '%f');
   end
end

%--------------------------------------------------------------------------

function t = create_table(ncol)
   t = struct('key' , [], ...
              'pos' ,  1, ...
              'data', zeros([0, ncol]));
end

%--------------------------------------------------------------------------

function [t, done] = append_data(t, data, ncol)
   nd   = numel(data);
   done = nd == 0;

   if ~done
      assert (mod(nd, ncol) == 1);

      cpty = size(t.data, 1);
      nrow = floor(nd / ncol);
      p    = t.pos(end) + nrow;

      if p > cpty
         % Extend storage space.
         t.data(2*p, 1 : ncol) = 0;
      end

      t.data(t.pos(end) : p-1, :) = reshape(data(2 : end), ncol, []) .';
      t.pos                       = [t.pos; p      ];
      t.key                       = [t.key; data(1)];
   end
end

%--------------------------------------------------------------------------

function t = finish_table(t)
   if numel(t.pos) == 1
      % Empty table.
      t.pos = [];
   else
      assert (t.pos(end) - t.pos(end-1) > 1,                 ...
             ['Table does not provide undersaturated data ', ...
              'in last record.']);

      % Trim any excess capacity.
      t.data = t.data(1 : t.pos(end) - 1, :);
   end
end

%--------------------------------------------------------------------------


function T = readRelPermTable(fid, kw, ntab, ncol)
%Input ECLIPSE/FrontSim saturation function table (SWOF/SGOF &c).
%
% SYNOPSIS:
%   Table = readRelPermTable(fid, kw, ntab, ncol)
%
% PARAMETERS:
%   fid  - Valid file identifier (as obtained using FOPEN) to file
%          containing the rel-perm table.  FTELL(fid) is assumed to be at
%          the start of the rel-perm table (i.e., after the keyword, if
%          any, which prompted the reading of this table).
%
%   kw   - Keyword that prompted this saturation function table read.  Used
%          in diagnostic messages only.
%
%   ntab - Number of relative permeability tables to input.  Typically
%          corresponds to the ECLIPSE/FrontSim parameter 'NTSFUN'.
%
%   ncol - Number of columns in each of the 'ntab' tables.
%
% RETURNS:
%   Table - An ntab-by-1 cell array of relative permeability tables.  Any
%           default values in columns 2:ncol specified by the syntax '1*'
%           in the input data are replaced by linearly interpolated data
%           from the non-default values in the corresponding column.
%
% SEE ALSO:
%   `swof`, `sgof`, `swfn`, `sgfn`.

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

   for t = 1 : ntab,
      % A rel-perm table is textually represented as a sequence of rows,
      % each row containing 'ncol' numbers (or defaulted elements of the
      % form '1*').  The sequence is terminated by a single slash ('/')
      % character.
      %
      % Read data into a character string, split on whitespace, convert to
      % DOUBLE, and (if required) insert default values computed by linear
      % interpolation.
      %
      data = readRecordString(fid);    % Treat table data as single record.
      data = splitString(removeQuotes(data));

      if isempty(data) || all(cellfun(@isempty, data))
         if t > 1,
            % Quoth Scripture:
            %
            %   The entire table may be defaulted provided the table is
            %   not the first.  Defaulted tables are replaced with a copy
            %   of the previous table.

            T{t} = T{t - 1};
         else
            error(msgid('FirstTable:Defaulted'), ...
                 ['The first table in saturation function ', ...
                  'keyword ''%s'' cannot be defaulted.'], kw);
         end
      else
         T{t} = convertTable(reshape(data, ncol, []) .');
         T{t} = insertDefaultValues(T{t}, t, ntab, kw);
      end
   end
end

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function T = convertTable(T)
   T = cellfun(@(s) sscanf(s, '%f'), stringReplace(T, '1*', 'NaN'));
end

%--------------------------------------------------------------------------

function T = insertDefaultValues(T, t, ntab, kw)
   for c = find(any(isnan(T))),
      assert (c ~= 1);
      i = isnan(T(:,c));

      if all(i),
         % Interpret an all-defaulted column as zero.  Happens, usually,
         % only for capillary pressure descriptions.  Warn about condition
         % to alert operator that this interpretation is in effect.
         %
         if ntab > 1,
            ndigits = floor(log10(ntab)) + 1;
            tabdsgn = sprintf(' in table %0*d/%0*d', ...
                              ndigits, t, ndigits, ntab);
         else
            tabdsgn = '';
         end

         warning(msgid('ColumnAllDefaulted'), ...
                ['All-defaulted column %d%s of saturation function ', ...
                 '''%s'' treated as all-zero.'], c, tabdsgn, kw);

         T(:,c) = 0;
      else
         T(i,c) = interp1(T(~i,1), T(~i,c), T(i,1), 'linear', 'extrap');
      end
   end
end

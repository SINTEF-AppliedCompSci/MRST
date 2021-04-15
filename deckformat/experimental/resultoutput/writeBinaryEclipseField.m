function fid = writeBinaryEclipseField(ecl_field, header, fid)
% Write an Eclipse field to a file using binary format.
%
% SYNOPSIS:
%   writeBinaryEclipseField(ecl_field, fid)
%
% PARAMETERS:
%   ecl_field - A structure with field data.
%   header    - A structure with the field''s name, size and data type.
%   fid       - Valid file identifier (as obtained using FOPEN) to file.
%
% RETURNS:
%   fid - File identifier.
%
% NOTE:
%   This is a fairly low-level function which should generally not be
%   invoked directly from user code.
%
% SEE ALSO:
%   `writeBinaryEclipseFile`.

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


  % Write heading
  fid = writeHeader(header, fid);

  % Write data array
  if(strcmp(header.data_type, 'CHAR'))
    fid = writeCharArray(ecl_field, header, fid);
  else
    fid = writeDataArray(ecl_field, header, fid);
  end
end


%%% =======================================================================
%   Private functions follow below
%
%%% -----------------------------------------------------------------------
function fid = writeHeader(header, fid)
  mfmt = 'ieee-be'; % Big endian machine format
  sep = 16;
  name = header.name;
  name(end+1:8) = ' ';
  fwrite(fid, sep,                 'int32', mfmt);  % start separator
  fwrite(fid, name,                'uint8', mfmt);  % keyword
  fwrite(fid, header.number,    'int32', mfmt);  % number of elements
  fwrite(fid, header.data_type, 'uint8', mfmt);  % data type (INTE,REAL,..
  fwrite(fid, sep,                 'int32', mfmt);  % end separator
end

%%% -----------------------------------------------------------------------
function fid = writeDataArray(ecl_field, header, fid)
  if (header.number == 0)
    return;
  end

  mfmt = 'ieee-be'; % Big endian machine format
  sep = 16;
  warg = setWarg(header.data_type);

  n = header.number;  % Number of elements
  max_block_size = warg.max_block_size;
  prec           = warg.precision;

  i=1;
  while (n > max_block_size)
    fwrite(fid, sep, 'int32', mfmt); % start separator
    ix = i:i+max_block_size-1;
    fwrite(fid, ecl_field.values(ix), prec, mfmt);
    fwrite(fid, sep, 'int32', mfmt); % end separator
    n = n - max_block_size;
    i = i + max_block_size;
  end
  fwrite(fid, sep, 'int32', mfmt);   % start separator
  fwrite(fid, ecl_field.values(i:i+n-1), prec, mfmt);
  fwrite(fid, sep, 'int32', mfmt);   % end separator
end

%%% -----------------------------------------------------------------------
function fid = writeCharArray(ecl_field, header, fid)
  if (header.number == 0)
    return;
  end

  mfmt = 'ieee-be'; % Big endian machine format
  sep = 16;
  warg = setWarg(header.data_type);

  % Convert cell array of strings to one character array.
   ecl_field.values{1}(end+1:8) = ' ';
   char_arr = char(ecl_field.values);
   ecl_field.values = reshape(char_arr',1, numel(char_arr));

  n = header.number;  % Number of elements
  max_block_size = warg.max_block_size;
  prec           = warg.precision;
  ecl_field.values = ecl_field.values';
  i=1;
  while (n > max_block_size)
    fwrite(fid, sep, 'int32', mfmt); % start separator
    ix = i:i+max_block_size*8-1;
    fwrite(fid, ecl_field.values(ix,:), prec, mfmt);
    fwrite(fid, sep, 'int32', mfmt); % end separator
    n = n - max_block_size;
    i = i + max_block_size*8;
  end
  fwrite(fid, sep, 'int32', mfmt);   % start separator
  fwrite(fid, ecl_field.values(i:i+n*8-1), prec, mfmt);
  fwrite(fid, sep, 'int32', mfmt);   % end separator
end

%%% -----------------------------------------------------------------------
function warg = setWarg(field_type)
switch field_type
  case 'INTE'
    warg.max_block_size = 1000;
    warg.precision = 'int32';
  case 'REAL'
    warg.max_block_size = 1000;
    warg.precision = 'float32';
  case 'DOUB'
    warg.max_block_size = 1000;
    warg.precision = 'float64';
  case 'CHAR'
    warg.max_block_size = 105;
    warg.precision = 'char';
  case 'LOGI'
    warg.max_block_size = 1000;
    warg.precision = 'int32';
  otherwise   % ???
    warg.max_block_size = 1000;
    warg.precision = 'int32';

end
end

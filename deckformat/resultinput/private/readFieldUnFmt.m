function [name, data] = readFieldUnFmt(fid, varargin)
%Extract next field from an unformatted ECLIPSE output file
%
% SYNOPSIS:
%   [name, output] = readFieldUnFmt(fid)
%
% PARAMETERS:
%   fid - Valid file identifier as obtained from FOPEN.  The file pointer
%         (FTELL(fid)) is assumed to be positioned directly before the next
%         keyword/field data.
%
% RETURNS:
%   name    - Name of this field.
%
%   output  - Field data.  A structure containing the following fields:
%              - type -- Data type of this field.
%                        String.  Supported values are
%                          - INTE -- Data values are Fortran INTEGERs.
%                          - REAL -- Data values are Fortran REALs.
%                          - DOUB -- Data values are Fortran DOUBLE PRECs.
%                          - LOGI -- Data values are Fortran LOGICALs.
%                          - CHAR -- Data values are character strings.
%
%              - values --
%                        Data values of this field.  A DOUBLE array for
%                        field data of type INTE, REAL, or DOUB.  A LOGICAL
%                        array for field data of type LOGI and a cell
%                        array of strings for field data of type CHAR.
%
% NOTE:
%   This is a fairly low-level function which should generally not be
%   invoked directly from user code.
%
% SEE ALSO:
%   `private/readFieldFmt`.

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

   assert (~isempty(fopen(fid)), ...
           'File identifier must point to open stream.');

   header = read_meta(fid);
   name   = header.name;
   ttype  = header.type;

   if header.ok
      switch lower(ttype)
         case {'inte', 'doub', 'real', 'logi'}
            values = read_record(fid, header);

            if any(strcmpi(ttype, {'DOUB', 'REAL'}))
               ttype = 'REAL';
            end

         case 'char'
            values = read_char_record(fid, header);

         otherwise
            if header.number == 0
               data = struct('type', ttype, 'values', {{}});
               return
            else
               error(msgid('Type:Unknown'), ...
                    ['Variable type ''%s'' is unexpected ', ...
                     'at this time.'], ttype);
            end
      end
   else
      values = {};
   end

   data = struct('type', ttype, 'values', values);
end

%--------------------------------------------------------------------------
% Private helpers follow.
%--------------------------------------------------------------------------

function header = read_meta(fid)

   name = fread(fid, 8, 'uint8=>char') .';

   eof = feof(fid);
   err = ferror(fid);

   if eof || ~isempty(err)
      name       = 'Input failure: ';

      if eof
         name = [name, 'End-of-file'];
      else
         name = [name, err];
      end

      nel        = -1;
      type       = 'FAIL';
      type_size  = -1;
      block_size = -1;
      type_prec  = NaN;
      ok         = false;

   else
      nel  = fread(fid, 1, 'int32', 'ieee-be');
      type = fread(fid, 4, 'uint8=>char') .';

      name = strtrim(name);

      if nel > 0
         [type_size, type] = get_type_characteristics(type, name);
         block_size        = get_block_size(type);
         type_prec         = get_type_prec (type, type_size, block_size);
      else
         type_size  = -1;
         type       = strtrim(type);
         block_size = -1;
         type_prec  = NaN;
      end

      ok = true;

      % Skip to start of record
      fseek(fid, 8, 'cof');
   end

   header = struct('number', nel       , ...
                   'type'  , type      , ...
                   'name'  , name      , ...
                   'size'  , type_size , ...
                   'blksiz', block_size, ...
                   'prec'  , type_prec , ...
                   'ok'    , ok        );
end

%--------------------------------------------------------------------------

function values = read_record(fid, header)
   if header.number > 0
      n      = header.number;
      prec   = header.prec;
      values = fread(fid, n, prec, 8, 'ieee-be');

      % Skip to start of next header unless n is a muliple of 1000 in which
      % case the skip has already been performed.
      if ~feof(fid) && (mod(n, 1000) ~= 0)
         fseek(fid, 8, 'cof');
      end

   else
      values = [];
   end
end

%--------------------------------------------------------------------------

function values = read_char_record(fid, header)
   if header.number > 0
      n      = header.number;
      nchar  = header.size;
      prec   = header.prec;
      skip   = 8;
      a      = reshape(fread(fid, n * nchar, prec, skip), 1, []);
      values = { cellstr(char(reshape(a, nchar, []) .')) };

      % Skip to start of next header unless 'n' is muliple of maximum block
      % size.  In the latter case, the FREAD call has already skipped over
      % the control characters.
      if ~feof(fid) && (mod(n, header.blksiz) ~= 0)
         fseek(fid, skip, 'cof');
      end
   else
      values = { [] };
   end
end

%--------------------------------------------------------------------------

function [size, type] = get_type_characteristics(type, name)
   if any(strcmp(type, {'INTE', 'REAL', 'LOGI'}))
      size = 4;

   elseif any(strcmp(type, {'CHAR', 'DOUB'}))
      size = 8;

   elseif ~isempty(regexp(type, '^C\d+$', 'once'))
      size = sscanf(type, 'C%d', 1);
      type = 'CHAR';  % Treat C0nn as CHAR for remainder of processing.

   else
      error('ElmType:UnSupp', ...
           ['Vector Element Type ''%s'' is not Supported ', ...
            'for Vector ''%s''.'], type, name);
   end
end

%--------------------------------------------------------------------------

function max_size = get_block_size(type)
   max_size = 1000;

   if strcmp(type, 'CHAR')
      max_size = 105;
   end
end

%--------------------------------------------------------------------------

function precision = get_type_prec(type, type_size, block_size)
   switch type
      case 'INTE'
         precision = sprintf('%d*int32', block_size);

      case 'REAL'
         precision = sprintf('%d*float32', block_size);

      case 'DOUB'
         precision = sprintf('%d*float64', block_size);

      case 'CHAR'
         % Usually just 105 * 8 == 840 bytes in a block, but this allows
         % for extended character types like C0nn too.
         precision = sprintf('%d*uchar=>char', block_size * type_size);

      case 'LOGI'
         precision = sprintf('%d*int32', block_size);

      otherwise
         error('ElmType:UnSupp', 'Field type ''%s'' is unsupported.', type);
   end
end

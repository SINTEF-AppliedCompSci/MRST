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
%   private/readFieldFmt.

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

   if header.ok,
      switch lower(ttype),
         case {'inte', 'doub', 'real', 'logi'},
            values = read_record     (fid, header);

            if any(strcmpi(ttype, {'DOUB', 'REAL'})),
               ttype = 'REAL';
            end

         case 'char',
            values = read_char_record(fid, header);

         otherwise,
            if header.number == 0,
               data = struct('type', ttype, 'values', {});
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
   skip_separator(fid);

   name = fread(fid, 8, 'uint8=>char') .';

   eof = feof(fid);
   err = ferror(fid);

   if eof || ~isempty(err);
      name       = 'Input failure: ';

      if eof,
         name    = [name, 'End-of-file'];
      else
         name    = [name, err];
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

      skip_separator(fid);

      name       = strtrim       (name);

      if nel > 0,
         type_size  = get_type_size (type);
         type_prec  = get_type_prec (type);
         block_size = get_block_size(type);
      else
         type_size  = -1;
         type       = strtrim(type);
         block_size = -1;
         type_prec  = NaN;
      end
      ok         = true;
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
   if header.number > 0,
      n      = header.number;
      prec   = header.prec;
      blksiz = header.blksiz;

      values = zeros([n, 1]);
      off    = 0;
      ix     = (1 : blksiz) .';

      while n > blksiz,
         values(off + ix) = readField(fid, blksiz, prec);

         n   = n   - blksiz;
         off = off + blksiz;
      end

      assert (off + n == numel(values), ...
              'Internal error in record reading.');

      values(off + (1 : n)) = readField(fid, n, prec);
   else
      values = [];
   end
end

%--------------------------------------------------------------------------

function values = read_char_record(fid, header)
   if header.number > 0,
      n      = header.number;
      prec   = header.prec;
      blksiz = header.blksiz;

      kws   = [];
      nchar = blksiz * header.size;
      while n > blksiz,
         a   = reshape(readField(fid, nchar, prec), 1, []);

         kws = [kws, a];         %#ok    % Willfully ignore MLINT warnings.

         n   = n - blksiz;
      end

      a      = reshape(readField(fid, n * header.size, prec), 1, []);
      values = { cellstr(char(reshape([kws, a], 8, []) .')) };
   else
      values = { [] };
   end
end

%--------------------------------------------------------------------------

function size = get_type_size(type)
   size = 4;

   if any(strcmp(type, {'CHAR', 'DOUB'})),
      size = 8;
   end
end

%--------------------------------------------------------------------------

function max_size = get_block_size(type)
   max_size = 1000;

   if strcmp(type, 'CHAR'),
      max_size = 105;
   end
end

%--------------------------------------------------------------------------

function precision = get_type_prec(type)
   switch type,
      case 'INTE', precision = 'int32';
      case 'REAL', precision = 'float32';
      case 'DOUB', precision = 'float64';
      case 'CHAR', precision = 'char';
      case 'LOGI', precision = 'int32';
      otherwise,
         error('Field type ''%s'' is unsupported.', type);
   end
end

%--------------------------------------------------------------------------

function values = readField(fid, n, precision)
   skip_separator(fid);

   values = fread(fid, n, precision, 'ieee-be');

   skip_separator(fid);
end

%--------------------------------------------------------------------------

function skip_separator(fid)
   fseek(fid, 4, 'cof');
end

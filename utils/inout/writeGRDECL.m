function varargout = writeGRDECL(grdecl, filename)
%Write a GRDECL structure out to permanent file on disk.
%
% SYNOPSIS:
%   writeGRDECL(grdecl, file)
%
% PARAMETERS:
%   grdecl - A corner-point GRDECL structure which may be passed to
%            function 'processGRDECL' in order to construct a grid
%            structure.
%
%   file   - Name of output file.  String.  This name is passed directly to
%            function 'fopen' using mode 'wt'.
%
% RETURNS:
%   Nothing.
%
% SEE ALSO:
%   `fopen`, `readGRDECL`, `processGRDECL`.

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

   success = true;

   [fid, msg] = fopen(filename, 'wt');
   if fid < 0, error(msg); end

   fprintf(fid, ['-- The following file was written from MATLAB ', ...
                 'using the MRST toolbox.\n\n']);

   if isfield(grdecl, 'cartDims'),
      fprintf(fid, 'SPECGRID\n%d %d %d 1 F\n/\n', grdecl.cartDims);
   end

   format = @(elmfmt, n) [repmat([' ', elmfmt], [1, n]), '\n'];

   write = @(elms_per_line, elmfmt, flds) ...
      write_fields(fid, grdecl, flds, format(elmfmt, elms_per_line));

   write( 6, '%12f'  , { 'COORD' });
   write( 8, '%12f'  , { 'ZCORN' });
   write(20, '%d'    , { 'ACTNUM', 'SATNUM' });
   write( 4, '%18.5e', [ { 'PORO' }, ...
                         strcat('PERM', { 'X', 'Y', 'Z' }) ]);

   fclose(fid);

   if nargout > 0,
      varargout{1} = success;
   end
end

%--------------------------------------------------------------------------

function write_fields(fid, grdecl, fields, fmt)
   for fld = reshape(fields(isfield(grdecl, fields)), 1, []),
      fprintf(fid, '\n%s\n', fld{1});
      fprintf(fid, fmt, grdecl.(fld{1}));
      fprintf(fid, '/\n');
   end
end

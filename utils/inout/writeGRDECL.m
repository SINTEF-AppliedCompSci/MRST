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
%   fopen, readGRDECL, processGRDECL.

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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

   format = @(s,n) [repmat([s, ' '], [1, n]), '\n'];

   [fid, msg] = fopen(filename, 'wt');
   if fid < 0, error(msg); end

   fprintf(fid, ['-- The following file was written from MATLAB ', ...
                 'using the MRST toolbox.\n\n']);

   if isfield(grdecl, 'cartDims'),
      fprintf(fid, 'SPECGRID\n%d %d %d 1 F\n/\n\n', grdecl.cartDims);
   end

   if isfield(grdecl, 'COORD')
      fprintf(fid, 'COORD\n');
      fprintf(fid, format('%12f', 6), grdecl.COORD);
      fprintf(fid, '/\n\n');
   end

   if isfield(grdecl, 'ZCORN'),
      fprintf(fid, 'ZCORN\n');
      fprintf(fid, format('%12f', 8), grdecl.ZCORN);
      fprintf(fid, '/\n\n');
   end

   if isfield(grdecl, 'ACTNUM'),
      fprintf(fid, 'ACTNUM\n');
      fprintf(fid, format('%d', 20), grdecl.ACTNUM);
      fprintf(fid, '/');
   end

   if isfield(grdecl, 'SATNUM'),
      fprintf(fid, '\n\nSATNUM\n');
      fprintf(fid, format('%d', 20), grdecl.SATNUM);
      fprintf(fid, '/');
   end

   if isfield(grdecl, 'PORO'),
      fprintf(fid, '\n\nPORO\n');
      fprintf(fid, format('%18.5e', 4), grdecl.PORO);
      fprintf(fid, '/');
   end

   if isfield(grdecl, 'PERMX'),
      fprintf(fid, '\n\nPERMX\n');
      fprintf(fid, format('%18.5e', 4), grdecl.PERMX);
      fprintf(fid, '/');
   end

   if isfield(grdecl, 'PERMY'),
      fprintf(fid, '\n\nPERMY\n');
      fprintf(fid, format('%18.5e', 4), grdecl.PERMY);
      fprintf(fid, '/');
   end

   if isfield(grdecl, 'PERMZ'),
      fprintf(fid, '\n\nPERMZ\n');
      fprintf(fid, format('%18.5e', 4), grdecl.PERMZ);
      fprintf(fid, '/');
   end

   fprintf(fid, '\n');
   fclose(fid);

   if nargout > 0,
      varargout{1} = success;
   end
end

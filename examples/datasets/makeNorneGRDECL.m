function ok = makeNorneGRDECL
%Create containing datafile for subset of Norne simulation model
%
% SYNOPSIS:
%   ok = makeNorneGRDECL
%
% DESCRIPTION:
%   This function ensures existence of the file
%
%       fullfile(getDatasetPath('norne'), 'NORNE.GRDECL')
%
%   that contains INCLUDE statements for an appropriate subset of the Norne
%   simulation model.
%
% PARAMETERS:
%   None.  This function operates on constant data.
%
% RETURNS:
%   ok - Whether or not we were able to create the datafile.
%
% SEE ALSO:
%   getDatasetPath.

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

   ok = grdecl_available();

   if (~ ok) && makeNorneSubsetAvailable(),
      ok = create_grdecl();
   end
end

%--------------------------------------------------------------------------

function ok = grdecl_available()
   ok = exist(grdecl_filename(), 'file') == 2;
end

%--------------------------------------------------------------------------

function ok = create_grdecl()
   assert (~ grdecl_available(), 'Internal Error');

   [fid, msg] = fopen(grdecl_filename(), 'wt');

   if fid < 0,
      warning('Open:Fail', 'Failed to open GRDECL output file: %s\n', msg);
   end

   write_grdecl(fid);

   fclose(fid);

   ok = true;
end

%--------------------------------------------------------------------------

function file = grdecl_filename()
   ndir = getDatasetPath('norne', 'skipAvailableCheck', true);
   file = fullfile(ndir, 'NORNE.GRDECL');
end

%--------------------------------------------------------------------------

function write_grdecl(fid)
   fprintf(fid, [ ...
'INCLUDE\n', ...
'  ''INCLUDE/IRAP_1005.GRDECL'' /\n\n', ...
'INCLUDE\n', ...
'  ''INCLUDE/ACTNUM_0704.prop'' /\n\n', ...
'INCLUDE\n', ...
'  ''INCLUDE/PERM_0704.prop'' /\n\n', ...
'INCLUDE\n', ...
'  ''INCLUDE/PORO_0704.prop'' /\n\n', ...
'INCLUDE\n', ...
'  ''INCLUDE/NTG_0704.prop'' /\n']);
end

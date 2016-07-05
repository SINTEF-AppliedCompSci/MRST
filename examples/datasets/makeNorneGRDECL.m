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
#COPYRIGHT#
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

function ok = makeNorneSubsetAvailable
%Ensure availability of subset of Norne simulation model
%
% SYNOPSIS:
%   ok = makeNorneSubsetAvailable
%
% DESCRIPTION:
%   If the dataset is not already present on disk in the directory named by
%
%       getDatasetPath('norne')
%
%   then this function will attempt to download a subset of the simulation
%   model from the OPM project's public dataset collection on GitHub
%   (https://github.com/OPM/opm-data/tree/master/norne).
%
% PARAMETERS:
%   None.  This function operates on constant data.
%
% RETURNS:
%   ok - Whether or not the Norne simulation model subset consisting of the
%        geometry description (corner-point format) and the petrophysical
%        properties (permeability, porosity, net-to-gross) are availble in
%        the Norne dataset directory.
%
% SEE ALSO:
%   getDatasetPath, githubDownload.

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

   [repo, rev, coll, files] = define_subset();

   idir      = include_directory();
   dfiles    = destination_files(idir, files);
   available = @() all(cellfun(@(x) exist(x, 'file') == 2, dfiles));

   ok = available();

   if ~ ok,
      ok = download(repo, rev, coll, files, idir) && available();
   end
end

%--------------------------------------------------------------------------

function [repo, rev, coll, files] = define_subset
   repo = 'OPM/opm-data';
   rev  = '2198d5b';
   coll = 'norne/INCLUDE';

   files = [ strcat('GRID/', { 'IRAP_1005.GRDECL' , ...
                               'ACTNUM_0704.prop' }), ...
             ...
             strcat('PETRO/', { 'PERM_0704.prop', ...
                                'PORO_0704.prop', ...
                                'NTG_0704.prop' }) ];
end

%--------------------------------------------------------------------------

function idir = include_directory
   ndir = getDatasetPath('norne', 'skipAvailableCheck', true);
   idir = fullfile(ndir, 'INCLUDE');
end

%--------------------------------------------------------------------------

function dfiles = destination_files(idir, files)
   [f, f, e] = cellfun(@fileparts, files, ...
                       'UniformOutput', false);                 %#ok<ASGLU>

   dfiles = fullfile(idir, strcat(f, e));
end

%--------------------------------------------------------------------------

function ok = download(repo, rev, coll, files, idir)
   sfiles = githubDownload(repo, 'Revision', rev, ...
                           'File', strcat(coll, '/', files));

   [ok, msg, id] = ensure_idir_exists(idir);

   if ok,
      ok = all(cellfun(@(f) movefile(f, idir), sfiles));
   else
      error(id, ['Unable to guarantee existence of Norne ', ...
                'INCLUDE directory: %s'], msg);
   end
end

%--------------------------------------------------------------------------

function [ok, msg, id] = ensure_idir_exists(idir)
   [ok, msg, id] = create_if_not_exists(idir);

   if ~ ok, return, end

   [ok, attr, id] = fileattrib(idir);

   if ~ ok,
      msg = sprintf(['Failed to obtain meta-data about ', ...
                     'directory (%s)'], attr);

   elseif ~ attr.UserWrite,
      id  = 'Directory:NotWritable';
      msg = 'Directory not user writable';
   end
end

%--------------------------------------------------------------------------

function [ok, msg, id] = create_if_not_exists(idir)
   [ok, msg, id] = deal(true, '', '');

   if ~ isdir(idir),
      [ok, msg, id] = mkdir(idir);

      if ~ ok,
         msg = sprintf('Failed to create directory (%s)', msg);
         return
      end
   end
end

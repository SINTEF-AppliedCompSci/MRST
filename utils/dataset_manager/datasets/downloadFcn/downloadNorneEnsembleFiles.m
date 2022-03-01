function ok = downloadNorneEnsembleFiles
%Ensure availability of subset of Norne simulation model
%
% SYNOPSIS:
%   ok = downloadNorneEnsembleFiles
%
% DESCRIPTION:
%   If the dataset is not already present on disk in the directory named by
%
%       getDatasetPath('norne')
%
%   then this function will attempt to download a subset of the simulation
%   model from the OPM project's public dataset collection on GitHub
%   (https://github.com/OPM/opm-tests/tree/master/norne).
%
% PARAMETERS:
%   None.  This function operates on constant data.
%
% RETURNS:
%   ok - Whether or not the Norne simulation model subset consisting of the
%        geometry description (corner-point format) and the petrophysical
%        properties (permeability, porosity, net-to-gross) is availble in
%        the Norne dataset directory.
%
% SEE ALSO:
%   `getDatasetPath`, `githubDownload`.

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

   [repo, rev, files] = define_subset();

   idir      = include_directory();
   dfiles    = destination_files(idir, files);
   available = @() all(cellfun(@(x) exist(x, 'file') == 2, dfiles));

   ok = available();

   if ~ ok
      ok = download(repo, rev, idir) && available();
   end

   if ok
      ok = copy_m_files(idir);
   end
end

%--------------------------------------------------------------------------

function [repo, rev, files] = define_subset
   repo = 'rolfjl/Norne-Initial-Ensemble';
   rev  = '56ed2be';
                            
    files = {'ACTNUM_0704.prop',...
                'FAULTMULT_AUG-2006.INC',...
                'NTG_0704.prop',...
                'PERM_0704.prop',...
                'PORO_0704.prop',...
                'FastGaussian.m',...
                'NorneGeostat.m'};      
end

%--------------------------------------------------------------------------

function idir = include_directory
   ndir = getDatasetPath('norne_ensemble', 'skipAvailableCheck', true);
   idir = fullfile(ndir, 'data');
end

%--------------------------------------------------------------------------

function dfiles = destination_files(idir, files)
   [f, f, e] = cellfun(@fileparts, files, ...
                       'UniformOutput', false);                 %#ok<ASGLU>

   dfiles = fullfile(idir, strcat(f, e));
end

%--------------------------------------------------------------------------

function ok = download(repo, rev, idir)
   sfiles = githubDownload(repo, 'Revision', rev, ...
                           'All',true,'Dest', idir);

   ok = ~ isempty(sfiles);
end

%--------------------------------------------------------------------------

function ok = copy_m_files(idir)
   scriptdir = fullfile(mrstPath('test-suite'), ...
                        'dataset-functions', 'norne');

   mfiles = { 'NorneGeostat.m', 'FastGaussian.m' };

   ok = true;
   for mfile = reshape(mfiles, 1, [])
      [stat, msg, id] = copyfile(fullfile(idir, mfile{1}), scriptdir);

      if stat ~= 1
         warning(id, ['Failed to copy supporting ensemble ', ...
                      'M-File ''%s'': %s'], mfile{1}, msg);

         ok = false;
      end
   end
end

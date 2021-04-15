function ok = makeSPE9OutputAvailable
%Ensure availability of the SPE 9 data set
%
% SYNOPSIS:
%   ok = makeSPE9OutputAvailable
%
% DESCRIPTION:
%   If the dataset is not already present on disk in the directory named by
%
%       getDatasetPath('spe9')
%
%   then this function will attempt to download a subset of the simulation
%   model from the OPM project's public dataset collection on GitHub
%   (https://github.com/OPM/opm-data/tree/master/spe9/eclipse-simulation).
%
% PARAMETERS:
%   None.  This function operates on constant data.
%
% RETURNS:
%   ok - Whether or not the requisite SPE-9 simulation output is availble
%        in the SPE-9 dataset directory.
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

   [repo, rev, coll, files] = define_subset();

   odir      = output_directory();
   dfiles    = destination_files(odir, files);
   available = @() all(cellfun(@(x) exist(x, 'file') == 2, dfiles));

   ok = available();

   if ~ ok
      ok = download(repo, rev, coll, files, odir) && available();
   end
end

%--------------------------------------------------------------------------

function [repo, rev, coll, files] = define_subset
   repo = 'OPM/opm-tests';
   rev  = '6c406ed';
   coll = 'spe9/eclipse-simulation';

   files = strcat('SPE9_CP', { '.INIT'   , ...
                               '.EGRID'  , ...
                               '.RSSPEC' , ...
                               '.SMSPEC' , ...
                               '.UNRST'  , ...
                               '.UNSMRY' });
end

%--------------------------------------------------------------------------

function odir = output_directory
   ndir = getDatasetPath('spe9', 'download', true, ...
                         'askBeforeDownload', false);
   odir = fullfile(ndir, 'Simulation-Output');
end

%--------------------------------------------------------------------------

function dfiles = destination_files(odir, files)
   [f, f, e] = cellfun(@fileparts, files, ...
                       'UniformOutput', false);                 %#ok<ASGLU>

   dfiles = fullfile(odir, strcat(f, e));
end

%--------------------------------------------------------------------------

function ok = download(repo, rev, coll, files, odir)
   bfiles = githubDownload(repo, 'Revision', rev, ...
                           'File', strcat(coll, '/', files), ...
                           'Dest', odir);

   ok = ~ isempty(bfiles);
end

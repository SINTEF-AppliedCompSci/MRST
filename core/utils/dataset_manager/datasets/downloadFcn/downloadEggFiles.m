function ok = downloadEggFiles()
%Ensure availability of the Egg Model dataset
%
% SYNOPSIS:
%   ok = downloadEggFiles
%
% DESCRIPTION:
%   If the dataset is not already present on disk in the directory named by
%
%       getDatasetPath('egg')
%
%   then this function will download the Egg Model archive from the SINTEF
%   website and extract the contents of the archive's top-level
%   'Egg_Model_Data_Files_v2' folder directly into that directory, rather
%   than leaving the files nested inside an extra 'Egg_Model_Data_Files_v2'
%   subfolder.
%
% PARAMETERS:
%   None.  This function operates on constant data.
%
% RETURNS:
%   ok - Whether or not the Egg Model dataset is available in the Egg
%        dataset directory.
%
% SEE ALSO:
%   `getDatasetPath`.

%{
Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.

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

   pth = getDatasetPath('egg', 'skipAvailableCheck', true);

   ok = dataset_present(pth);

   if ~ok
      ok = download_and_extract(pth) && dataset_present(pth);
   end
end

%--------------------------------------------------------------------------

function tf = dataset_present(pth)
   tf = isdir(pth) && numel(dir(pth)) > 2; %#ok<ISDIR>  (ignore '.' and '..')
end

%--------------------------------------------------------------------------

function ok = download_and_extract(pth)
   url = 'https://www.sintef.no/globalassets/project/mrst/egg.zip';

   if ~isdir(pth)
      mkdir(pth);
   end

   tmpdir = tempname();
   mkdir(tmpdir);

   ok = false;
   try
      wget    = mrstWebSave();
      zipfile = fullfile(tmpdir, 'egg.zip');

      if ~isempty(wget(zipfile, url))
         unzip(zipfile, tmpdir);

         srcdir = fullfile(tmpdir, 'Egg_Model_Data_Files_v2');

         if isdir(srcdir)
            copy_contents(srcdir, pth);

            ok = true;
         end
      end
   catch me
      warning('downloadEggFiles:Failed', ...
               'Failed to download Egg dataset: %s', me.message);
   end

   if isdir(tmpdir)
      rmdir(tmpdir, 's');
   end
end

%--------------------------------------------------------------------------

function copy_contents(srcdir, pth)
   items = dir(srcdir);
   items = items(~ismember({ items.name }, { '.', '..' }));

   for i = 1 : numel(items)
      [status, msg, id] = copyfile(fullfile(srcdir, items(i).name), ...
                                    fullfile(pth,    items(i).name));

      if status ~= 1
         warning(id, 'Failed to copy ''%s'': %s', items(i).name, msg);
      end
   end
end

function varargout = downloadDataset(name, askFirst)
%Download a dataset given by name, subject to availability
%
% SYNOPSIS:
%   [pth, ok] = downloadDataset('datasetname')
%   % Do not ask for download permission
%   [pth, ok] = downloadDataset('datasetname', false)
%
% DESCRIPTION:
%   Download a dataset, if available. If the dataset is not publicly
%   available due to technical or license related reasons, instructions
%   for manual download (if any) along with the dataset website will be
%   printed.
%
% REQUIRED PARAMETERS:
%   name     - The name of the dataset.
%
%   askFirst - (OPTIONAL) Will ask in the command window before downloading
%              any files, along with a estimate of how large the file(s)
%              to be downloaded are.
%
%
% RETURNS:
%   pth      - Path where dataset was / should be placed.
%
%   ok       - Boolean indicating if the download was sucessful.
%
% SEE ALSO:
%   getAvailableDatasets, mrstDatasetGUI

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

    if nargin == 1
        askFirst = true;
    end
    
    [info, present] = feval(['dataset_', lower(name)]);
    pth = fullfile(mrstDataDirectory(), info.name);
    
    if present
        dispif(mrstVerbose(), ...
               'Dataset already installed. Current path:\n %s\n', pth);
        ok = true;

    elseif datasetHasCustomDownloadFcn(info)
        % Dataset provides custom download function.  We assume that the
        % download function knows how to make the dataset available in
        % 'pth' so we'll just call it without any arguments.

        ok = do_download(info, askFirst, info.downloadFcn{1});

    elseif ~datasetHasValidFileURL(info)
        % Empty (or non-string) file URL.  Don't know how to do this.

        fail(info);
        ok = false;
    else
        % Static, non-empty string URL for dataset file.  Direct Download.

        ok = do_download(info, askFirst, @()download_URL(info, pth));
    end

    if nargout
        varargout{1} = pth;
        if nargout > 1
            varargout{2} = ok;
        end
    end
end

%--------------------------------------------------------------------------

function s = check_case_sensitive
   if ispc,
      s = 'case insensitive';
   else
      s = 'case sensitive';
   end
end

%--------------------------------------------------------------------------

function fail(info)
   % Display instructions
   fprintf('Dataset ''%s'' not directly available for download.\n', ...
           info.name);

   if ~isempty(info.instructions)
      fprintf(['*** Instructions for downloading:\n', ...
               '*** %s\n'], info.instructions);
   end

   if ~isempty(info.website)
      fprintf(['Go to the dataset webpage for more information:\n\n', ...
               '   <a href="%s">%s</a>\n\n'], info.website, info.website);
   end

   fprintf(['Dataset should be placed in directory:\n\n', ...
            '   %s (%s)\n'], pth, check_case_sensitive());
end

%--------------------------------------------------------------------------

function ok = download_URL(info, pth)
   furl = info.fileurl;

   assert (ischar(furl) && ~ isempty(furl), ...
           'Internal Error in Dataset Download From Static URL');

   [ext, ext, ext] = fileparts(furl);                           %#ok<ASGLU>

   fprintf('Attempting to download ''%s'' dataset (%.1f MB)...\n', ...
           info.name, info.filesize);

   is_known = any(strcmpi(ext, { '.gz', '.tgz', '.zip' }));

   if ~ is_known,
      warning('URLType:Unknown', ...
              'Unknown (URL) file extension: %s', ext);

      ok = false;
   else
      switch lower(ext)
         case { '.gz', '.tgz' },
            untar(furl, pth);

         case '.zip'
            unzip(furl, pth);

         otherwise
            error('Internal Error: Condition Mismatch');
      end

      ok = true;
      fprintf('Successfully downloaded dataset.\n');
   end
end

%--------------------------------------------------------------------------

function ok = do_download(info, askFirst, downloadFcn)
   do_get = true;
   if askFirst
      do_get = confirm_download(info);
   end

   ok = false;
   if do_get
      ok = downloadFcn();
   end
end

%--------------------------------------------------------------------------

function ok = confirm_download(info)
   prompt   = sprintf(['Do you want to download dataset ''%s'' ', ...
                       '(%.2f MB), Y/N [Y]:'], info.name, info.filesize);
   response = input(prompt, 's');

   ok = isempty(response) || any(strcmpi(response, {'y', 'yes'}));
end

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
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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
    
    [info, present] = eval(['dataset_', lower(name), '()']);
    name = info.name;
    pth = fullfile(mrstDataDirectory(), name);
    
    if present
        dispif(mrstVerbose(), 'Dataset already installed. Current path:\n %s\n', pth);
        ok = true;
    else
        furl = info.fileurl;
        if isempty(furl)
            % Display instructions
            disp(['Dataset ''', name, ''' not directly available for download.']);
            if ~isempty(info.instructions)
                fprintf('\n*** Instructions for downloading: \n*** %s \n', info.instructions);
            end
            if ~isempty(info.website)
                fprintf('\nGo to the dataset webpage for more information: \n\n    <a href="%s">%s</a> \n \n', info.website, info.website);
            end
            fprintf('Dataset should be placed in directory:\n\n    %s (%s)\n', pth, check_case_sensitive());
            ok = false;
        else
            % We know how to download it
            [tmp, tmp, ext] = fileparts(furl); %#ok
            if askFirst
                str = sprintf('Do you want to download dataset ''%s'' (%s MB), Y/N [Y]:', name, num2str(info.filesize));
                v = input(str,'s');
                switch(lower(v))
                    case {'y', '', 'yes'}
                        
                    otherwise
                        disp('Aborting. Dataset not downloaded.');
                        return
                end
            else
                
            end
            fprintf('Attempting to download ''%s'' dataset (%1.1f MB)...\n', ...
                                name, info.filesize);
            switch lower(ext)
                case { '.gz', '.tgz' },
                    untar(furl, pth);
                case '.zip'
                    unzip(furl, pth);
                otherwise
                    error(['Unknown file extension: ', ext]);
            end
            fprintf('Successfully downloaded dataset.\n');
            ok = true;
        end
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

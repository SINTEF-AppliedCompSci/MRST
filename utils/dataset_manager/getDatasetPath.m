function pth = getDatasetPath(name, varargin)
%Get the path of a dataset (optionally: try to download it if missing)
%
% SYNOPSIS:
%   pth = getDatasetPath('spe1');
%   pth = getDatasetPath('spe1', 'download', true);
%
% REQUIRED PARAMETERS:
%   name    -  The name of the dataset. Must be known to MRST, see
%              getAvailableDatasets for details.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   download - Boolean indicating if the functions should attempt to
%              download the dataset if it is missing. Default: Enabled.
%
%   askBeforeDownload - Boolean. If the download option is enabled, setting
%              this to true will prompt the user before starting a
%              potentially large download. Default: Enabled.
% RETURNS:
%   pth      - Path to dataset.
%
% NOTE: 
%   If this function returns, the dataset will be present at the path
%   given. Any other situation will result in an error being thrown.
% 
% SEE ALSO:
%   mrstdatasetGUI, downloadDataset

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
    opt = struct('download',          true, ...
                 'askBeforeDownload', true);
    opt = merge_options(opt, varargin{:});
    [info, present] = getDatasetInfo(name);

    if present
        % Return path to dataset
        pth = fullfile(mrstDataDirectory(), info.name);
    else
        if opt.download
            % We want to attempt a download
            [pth, ok] = downloadDataset(name, opt.askBeforeDownload);
            if ~ok
                error(['Unable to load dataset ''', info.name,...
                       ''' due to aborted/missing download']);
            end
        else
            % Nothing more to do, throw error
            error(['Dataset ''', info.name, ''' not found. ', ...
                   'Consider running ''downloadDataset(''', info.name, ''')''']);
        end
    end 
end
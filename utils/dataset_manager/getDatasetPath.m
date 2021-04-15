function pth = getDatasetPath(name, varargin)
%Get the path of a dataset (optionally: try to download it if missing)
%
% SYNOPSIS:
%   pth = getDatasetPath(name)
%   pth = getDatasetPath(name, 'pn1', pv1, ...)
%
% REQUIRED PARAMETERS:
%   name    -  The name of the dataset. Must be known to MRST, see
%              `getAvailableDatasets` for details.
%
% OPTIONAL PARAMETERS:
%   download - Boolean indicating if the functions should attempt to
%              download the dataset if it is missing. Default: Enabled.
%
%   askBeforeDownload - Boolean. If the download option is enabled, setting
%                       this to true will prompt the user before starting a
%                       potentially large download. Default: Enabled.
%
%   skipAvailableCheck -
%              Boolean flag indicating whether or not to omit checking for
%              presence of data on disk.  This is mainly intended for the
%              case of needing to manually download objects into the
%              dataset's containing directory through some external means.
%              Default value: `skipAvailableCheck = false` (*do* check if
%              the dataset is available).
%
% RETURNS:
%   pth      - Path to dataset.
%
% NOTE: 
%   If this function returns, the dataset will be present at the path
%   given. Any other situation will result in an error being thrown.
%
%   Using `skipAvailableCheck` bypasses this basic safety measure of MRST's
%   dataset handling.  Consequently, said option should be used only when
%   circumstances so dictate.
% 
% SEE ALSO:
%   `mrstdatasetGUI`, `downloadDataset`.

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

    opt = struct('download',           mrstSettings('get', 'allowDL'), ...
                 'askBeforeDownload',  mrstSettings('get', 'promptDL'), ...
                 'skipAvailableCheck', false);
    opt = merge_options(opt, varargin{:});
    [info, present] = getDatasetInfo(name);

    if present || opt.skipAvailableCheck
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

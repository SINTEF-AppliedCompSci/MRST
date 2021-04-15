function downloadAllDatasets(askFirst)
%Download all datasets known to MRST and available for direct download
%
% REQUIRED PARAMETERS:
%   
%   askFirst - (OPTIONAL) Boolean indicating if the user should be prompted
%              before starting the download. Default on. 
%
% RETURNS:
%   Nothing. Prints output to command line.
%
% NOTE:
%   Avoid using this function in scripts. Rather, opt for the more
%   conservative `downloadDataset` function to have control over which
%   files to download.
%
% SEE ALSO:
%   `mrstDatasetGUI`, `downloadDataset`

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


    if nargin == 0
        askFirst = true;
    end
    [info, present] = getAvailableDatasets();

    info = info(~present);

    pick = arrayfun(@datasetHasValidFileURL     , info) | ...
           arrayfun(@datasetHasCustomDownloadFcn, info);

    info = info(pick);
    
    if isempty(info)
        disp('All known datasets are already installed');
        return
    end
    if askFirst
        str = sprintf('Do you want to download ALL %d datasets (%d MB), Y/N [Y]:',...
                        numel(info), ceil(sum([info.filesize])));
        v = input(str,'s');
        switch(lower(v))
            case {'y', '', 'yes'}
                % Do nothing, user accepted
            otherwise
                disp('Aborting. No downloads started.');
                return
        end
    end
    
    % Download all datasets
    timer = tic();
    nI = numel(info);
    ok = false(nI, 1);
    for i = 1:nI
        [tmp, ok(i)] = downloadDataset(info(i).name, false); %#ok
    end
    T = toc(timer);
    % Note to the user how it went
    fprintf('Downloaded %d files in %s', sum(ok), formatTimeRange(T));
    if ~all(ok)
        fprintf('%d download(s) failed', sum(ok));
    end
    fprintf('\n');
end

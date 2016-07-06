function [info, present] = getAvailableDatasets()
%Get a list of structs indicating possible and present datasets in MRST
%
% SYNOPSIS:
%   [info, present] = getAvailableDatasets
%
% DESCRIPTION:
%   This function gathers info structs about all datasets MRST knows about.
%   The info struct will contain useful information about the dataset, such
%   as the name, type of model, download location, any MRST examples using
%   that dataset and so on.
%
%   getAvailableDatasets relies on functions in the "datasets" directory o
%   nthe form "dataset_name" which produces the info structs and checks if
%   the dataset is already present.
%
% REQUIRED PARAMETERS:
%   None.
%
% RETURNS:
%   info    - A array of structs, where each element corresponds to
%             information about a given dataset. See datasetInfoStruct for
%             possible fields, along with their meanings.
% 
%   present - Boolean array of the same size as info. If present(i) is
%             true, the dataset with name info(i).name is already
%             downloaded and available for use by MRST.
% SEE ALSO:
%   mrstDatasetGUI, datasetInfoStruct

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
    pth = mfilename('fullpath');
    pth = pth(1:(end-numel(mfilename()) - 1));
    sets = dir(fullfile(pth, 'datasets'));
    
    % Extract functions that start with dataset and is not directories
    sn = 'dataset';
    ok = strncmp({sets.name}, sn, numel(sn)) & ~[sets.isdir];

    names = {sets(ok).name}';
    
    nD = numel(names);
    
    % Loop over candidate functions, adding to array as we go (dynamic
    % expansion is ok since the number of datasets is relatively small).
    [info, present] = deal([], logical([]));
    for i = 1:nD
        [name, name, ext] = fileparts(names{i}); %#ok
        if ~strcmpi(ext, '.m')
            continue
        end
        
        [info_i, ok] = eval([name, '()']);
        info = [info; info_i];  %#ok
        present = [present; ok];%#ok
    end
end

function [info, present] = getDatasetInfo(name)
% Get info struct for a given dataset.
%
% SYNOPSIS:
%   [info, present] = getDatasetInfo('datasetname')
%
% REQUIRED PARAMETERS:
%   name    - Dataset name. Must be known to MRST.
%
%
% RETURNS:
%   info    - Info struct as defined by datasetInfoStruct with containing
%             metadata about the dataset with the supplied name.
%
% SEE ALSO:
%   `datasetInfoStruct`, `getDatasetPath`, `listDatasetExamples`

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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


    dname = ['dataset_', lower(name)];
    if exist(dname, 'file') == 2
        [info, present] = eval([dname, '()']);
    else
        error(['Dataset ''', name, ''' is not known to MRST.']);
    end
end
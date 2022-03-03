function [info, present] = getDatasetInfo(name)
%Retrieve Descriptive Information About a Named Dataset
%
% SYNOPSIS:
%   info            = getDatasetInfo(name)
%   [info, present] = getDatasetInfo(...)
%
% PARAMETERS:
%   name - Dataset name.  Character vector or, if supported in the running
%          MATLAB version, a string type.  Name must be known to MRST.
%
% RETURNS:
%   info    - Descriptive structure defined by function `datasetInfoStruct`
%             containing information about the named dataset.
%
%   present - Whether or not physical files pertaining to the named dataset
%             already exist on the local computer system's disk in MRST's
%             managed dataset location.
%
% SEE ALSO:
%   `datasetInfoStruct`, `getDatasetPath`, `listDatasetExamples`.

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

    dname = strcat('dataset_', lower(name));
    if any(exist(dname, 'file') == [2, 3, 6]) % M, MEX, or P in path
        [info, present] = feval(dname);
    else
        error(['Dataset ''', name, ''' is not known to MRST.']);
    end
end

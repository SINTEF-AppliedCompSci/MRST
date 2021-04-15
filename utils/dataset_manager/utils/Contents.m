% UTILS
%   Supporting routines for managing collection of external datasets
%
% Files
%   datasetInfoStruct           - Get a struct containing standardized information about a dataset
%   datasetHasCustomDownloadFcn - Predicate for whether or not a dataset provides a custom download function
%   datasetHasValidFileURL      - Predicate for whether or not a dataset provides a valid file URL
%   downloadAllDatasets         - Download all datasets known to MRST and available for direct download
%   getDatasetInfo              - Get info struct for a given dataset.
%   listDatasetExamples         - List all MRST examples using a specific dataset

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

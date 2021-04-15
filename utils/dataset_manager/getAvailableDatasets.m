function [info, present] = getAvailableDatasets()
%Get a list of structures indicating possible and present datasets in MRST
%
% SYNOPSIS:
%   [info, present] = getAvailableDatasets
%
% DESCRIPTION:
%   This function collects descriptive information about all datasets known
%   to MRST.  One structure contains useful information about a particular
%   dataset, such as its name, type of model (e.g., characteristics of the
%   model geometry, the model's petrophysical properties, active fluid
%   phases, common simulation scenarios), online location (URL) from which
%   to download the data, any MRST examples using that dataset and so on.
%
%   Function `getAvailableDatasets` relies on dedicated helper functions
%   named `dataset_*` (lower-case) located in the sudirectory `datasets` of
%   this function's location.  Each helper function produces a single
%   collection of metadata about a particular dataset and determines if the
%   associated dataset is already present (i.e., downloaded).
%
%   Creating a new dataset consequently means implementing a new
%   `dataset_*` function.
%
% PARAMETERS:
%   None.
%
% RETURNS:
%   info    - An array of structures in which each array element provides
%             information about a separate dataset.  See datasetInfoStruct
%             for possible fields, along with their meanings.
% 
%   present - Logical array of the same size as `info`.  If `present(i)` is
%             true, the dataset named `info(i).name` is already downloaded
%             and exists in the directory returned by function ::
%
%                 mrstDataDirectory
%
%             This dataset is therefore directly available for use by MRST.
%             Otherwise, the dataset may be retrieved interactively through
%             a graphical user interface (function `mrstDatasetGUI`) or
%             programmatically by means of function `downloadDataset` or
%             `downloadAllDatasets`.
%
% SEE ALSO:
%   `mrstDatasetGUI`, `mrstDataDirectory`, `datasetInfoStruct`, `downloadDataset`.

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

   % Loop over candidate functions, adding to arrays as we go (dynamic
   % expansion is ok since the number of datasets is relatively small).
   [info, present] = deal([], logical([]));

   for dataset = reshape(identify_datasets(), 1, []),
      [info, present] = inspect_dataset(info, present, dataset{1});
   end
end

%--------------------------------------------------------------------------

function sets = identify_datasets()
   sets = what(fullfile(this_directory(), 'datasets'));
   sets = sets.m(is_dataset(sets.m));

   edir = fullfile(this_directory(), 'datasets', 'experimental');
   if isdir(edir)
      s    = what(edir);
      sets = [ sets ; s.m(is_dataset(s.m)) ];
   end
end

%--------------------------------------------------------------------------

function [info, present] = inspect_dataset(info, present, dataset)
   [name, name] = fileparts(dataset);                           %#ok<ASGLU>

   [info_i, ok] = feval(name);

   info    = [info    ; info_i];
   present = [present ; ok    ];
end

%--------------------------------------------------------------------------

function d = this_directory()
   d = fileparts(mfilename('fullpath'));
end

%--------------------------------------------------------------------------

function tf = is_dataset(sets)
   tf = ~ cellfun(@isempty, regexp(sets, '^dataset_'));
end

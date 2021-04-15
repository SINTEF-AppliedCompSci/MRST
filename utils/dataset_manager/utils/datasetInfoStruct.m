function [info, present] = datasetInfoStruct(varargin)
%Get a struct containing standardized information about a dataset
%
% SYNOPSIS:
%   [info, present] = datasetInfoStruct('name', 'myDataset', ...)
%
% OPTIONAL PARAMETERS:
%
%   'name'          - Name of the dataset. Any capitalization will be kept,
%                     but MRST does not allow multiple datasets with the
%                     same name aside from capitalization. Any whitespace
%                     will be stripped from the name.
%   'description'   - Description of the dataset.
%   'cells'         - Number of cells in the dataset, if applicable.
%   'website'       - Website for the dataset. Should be the home page of
%                     the dataset, if it exists.
%   'fileurl'       - URL to directly download the dataset. Supports files
%                     in standard archive formats (.tar.gz, .tgz, .zip).
%   'filesize'      - Size of files (in MegaBytes). Used to give the user a
%                     estimate before downloading.
%   'instructions'  - If fileurl is not provided, it typically means that
%                     the dataset is not directly available. The
%                     instructions field should tell the user how to get
%                     the dataset in that case (or if it is not available
%                     at all).
%   'hasGrid'       - Boolean. True if the dataset includes a grid.
%   'hasRock'       - Boolean. True if the dataset includes petrophysical
%                     data for the rock. 
%   'hasFluid'      - Boolean. True if the dataset includes a fluid model.
%   'source'        - String indicating where the dataset originated.
%   'examples'      - Cell array with strings of the form
%                     `module:examplename`. For example, setting it to::
%
%                       {'module_a:test1',...
%                        'module_a:test5', ...
%                        'module_b: myTest'};
%
%                     will indicate that the examples in module_a named
%                     test1 and test5 use the dataset and that the test
%                     myTest in module_b also does use it.
%                     
%   'modelType'     - The type of the model (e.g. grid only, black oil,
%                     corner point).
%
%   'downloadFcn'   - Custom dataset download function.  Handle to function
%                     that takes no input parameters and returns a status
%                     code signifying whether or not the file data was
%                     successfully downloaded.  If specified, this value
%                     takes precedence over the `fileurl` parameter.
%
%   'note'          - Additional notes concerning the dataset, e.g.
%                     special actions to take before using the datafiles.
%                     Character vector or cell array of character vectors
%                     ('cellstring').  Default value: `note = ''` (no
%                     additional notes).
%
% RETURNS:
%   info            - Struct with info fields as listed above.
%
%   present         - Boolean indicator if the directory 
%                     `fullfile(mrstDataDirectory(), info.name)` exists.
% SEE ALSO:
%   `mrstDatasetGUI`

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

    info =   struct('name',                '', ...
                    'description',         '', ...
                    'cells',               [], ...
                    'website',             '', ...
                    'fileurl',             '', ...
                    'filesize',            [], ...
                    'instructions',        '', ...
                    'hasGrid',             false, ...
                    'hasRock',             false, ...
                    'hasFluid',            false, ...
                    'source',              '', ...
                    'image',               '', ...
                    'examples',            {{}}, ...
                    'modelType',           'Unknown', ...
                    'downloadFcn',         {{}}, ...
                    'note',                {{}} ...
                    );
    info = merge_options(info, varargin{:});
    info.name = info.name(~isspace(info.name));
    
    present = isdir(fullfile(mrstDataDirectory(), info.name));

    % Work out if a image is attached to the dataset
    pth = fileparts(mfilename('fullpath'));
    imgdir = fullfile(pth, '..', 'datasets', 'img');
    imgs = dir(imgdir);
    
    imgNames = {imgs.name};
    match = regexpi(imgNames,[lower(info.name), '.(png|jpg|jpeg|gif)'],'match');
    match = match(~cellfun(@isempty, match));
    if ~isempty(match)
        info.image = fullfile(imgdir, match{1}{1});
    end
end

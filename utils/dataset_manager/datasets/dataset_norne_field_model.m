function [info, present] = dataset_norne_field_model()
% Info function for Norne dataset. Use getDatasetInfo or getAvailableDatasets for practical purposes.

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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


    [info, present] = datasetInfoStruct(...
        'name', 'norne_field_model', ...
        'fileurl', 'https://www.sintef.no/globalassets/project/mrst/norne_field_model_example_files.zip', ...
        'hasGrid', true, ...
        'hasRock', true, ...
        'hasFluid', true, ...
        'filesize', 69.0, ...
        'examples', {'ad-blackoil:fieldModelNorneExample'}, ...
        'description', ['Datafiles required to run the fieldModelNorne example'], ...
        'modelType', 'Three-phase, black-oil, corner-point' ...
         );
end

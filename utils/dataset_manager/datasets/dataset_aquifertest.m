function [info, present] = dataset_aquifertest()
% Info function for the dataset of a 2D oil-water two phase cases with Fetkovich aquifers. Use getDatasetInfo or getAvailableDatasets for practical purposes.

%{
Copyright 2015-2021 SINTEF Digital, Mathematics & Cybernetics.

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
        'name', 'aquifertest', ...
        'fileurl', 'https://www.sintef.no/contentassets/124f261f170947a6bc51dd76aea66129/aquifertestdata.zip', ...
        'hasGrid', true, ...
        'hasRock', true, ...
        'hasFluid', true, ...
        'examples', {'ad-blackoil:aquifertest'}, ...
        'description', 'This is a simple model with aquifer to test aquifer functionality.', ...
        'modelType', 'Two-phase with aquifer. Cartesian grid' ...
         );
end

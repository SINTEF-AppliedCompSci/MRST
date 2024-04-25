function [info, present] = dataset_co2lab_sampled_tables()
% Info function for CO2lab sampled tables dataset. Use getDatasetInfo or getAvailableDatasets for practical purposes.

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
        'name', 'CO2lab_Sampled_Tables', ...
        'website', '', ...
        'fileurl', 'https://www.sintef.no/contentassets/124f261f170947a6bc51dd76aea66129/co2lab_sampled_tables.zip', ...
        'examples', {}, ...
        'description', 'A set of precomputed sampled tables for the CO2lab module.', ...
        'filesize',    102.8, ...
        'modelType', 'MAT files with tables' ...
         );
end

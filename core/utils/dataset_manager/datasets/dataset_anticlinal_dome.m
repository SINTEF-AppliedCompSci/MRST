function [info, present] = dataset_anticlinal_dome()
    % Info function for geothermal dataset. Use getDatasetInfo or getAvailableDatasets for practical purposes.
    
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
            'name'       , 'H2_STORAGE_RSRV_DOME', ...
            'description', 'Data for the h2store module examples `simulateStorageAnticlinalDomeAquifer`', ...
            'fileurl'    , 'https://www.sintef.no/contentassets/124f261f170947a6bc51dd76aea66129/anticlinal_dome.zip', ...
            'hasGrid'    , true, ...
            'hasRock'    , true, ...
            'hasFluid'   , true, ...
            'examples'   , {'simulateStorageAnticlinalDomeAquifer.m'}, ...
            'filesize'   , 0.159 ...
            );
    end
    
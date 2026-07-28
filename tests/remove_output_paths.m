%% Copyright notice
%{
Copyright 2026, NORCE Research AS, Computational Geosciences and
Modeling.

This file is part of the ad-micp module.

ad-micp is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ad-micp is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file. If not, see <http://www.gnu.org/licenses/>.
%}

function remove_output_paths(outputDirectory)
    pathEntries = strsplit(path(), pathsep());

    for entryIndex = 1 : numel(pathEntries)
        entry = pathEntries{entryIndex};

        isOutputDirectory = strcmp(entry, outputDirectory);
        isOutputSubdirectory = strncmp( ...
            entry, [outputDirectory filesep()], ...
            numel(outputDirectory) + 1);

        if (isOutputDirectory || isOutputSubdirectory) && ...
                isfolder(entry)
            rmpath(entry);
        end
    end
end

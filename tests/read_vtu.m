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

function data = read_vtu(filename)
    txt = fileread(filename);
    data = struct();
    % Find all DataArray names in the file
    names = regexp(txt, ...
        '<DataArray[^>]*Name="([^"]+)"', ...
        'tokens');
    names = unique(cellfun(@(x) x{1}, names, ...
        'UniformOutput', false));
    for i = 1:numel(names)
        name = names{i};
        pattern = sprintf( ...
            '<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>', ...
            regexptranslate("escape", name));
        token = regexp(txt, pattern, ...
                       "tokens", ...
                       "once", ...
                       "dotexceptnewline");
        if isempty(token)
            token = regexp(txt, pattern, ...
                           "tokens", ...
                           "once");
        end
        if isempty(token)
            continue;
        end
        data.(name) = sscanf(token{1}, "%f");
    end
end

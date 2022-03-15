function settings = struct2settings(oldsettings)
%Undocumented Utility Function

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    settings = settingsStruct();
    fnames = fieldnames(oldsettings);

    supported = isfield(settings, fnames);
    if ~all(supported)
       unsupported = sprintf(' * %s\n', fnames{~supported});
       pl = ''; if sum(~supported) ~= 1, pl = 's'; end
       error('Unsupported setting name%s\n%s\nCannot convert to ''settingsStruct''.', pl, unsupported);
    end

    for setting = reshape(fnames, 1, [])
       settings.(setting{1}) = oldsettings.(setting{1});
    end
end

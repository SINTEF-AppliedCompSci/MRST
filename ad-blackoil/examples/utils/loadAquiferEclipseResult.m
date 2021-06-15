function [statesEcl, G] = loadAquiferEclipseResult(dir, fnroot)
%Undocumented Utility Function

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

    fn = fullfile(dir, [fnroot, '.INIT']);
    init = readEclipseOutputFileUnFmt(fn);
    fn = fullfile(dir, [fnroot, '.EGRID']);
    grid = readEclipseOutputFileUnFmt(fn);
    [eclG, eclrock, N, T] = initGridFromEclipseOutput(init, grid, 'outputSimGrid', ...
                                                      true);
    G = eclG{1};

    statesEcl = {};
    for i = 1 : 21
        fnname = sprintf('%s.X%04.0f', fnroot, i);
        fn = fullfile(dir, fnname);
        stateEclFmt = readEclipseOutputFileUnFmt(fn);
        stateEcl.pressure = (stateEclFmt.PRESSURE.values);
        sw = stateEclFmt.SWAT.values;
        stateEcl.s = [sw, (1 - sw)];
        statesEcl{i} = stateEcl;
    end
end

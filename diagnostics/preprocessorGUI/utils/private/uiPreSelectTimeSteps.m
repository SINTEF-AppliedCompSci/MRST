function steps = uiPreSelectTimeSteps(info)
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

d = dialog('Resize', 'on');
s = PreDiagnosticsSelector('Parent', d, 'restartInfo', info);
d.CloseRequestFcn = @setSteps;
d.SizeChangedFcn  = @updateSize;
uiwait(d)

    function setSteps(src, event)
    steps = s.ix;
    delete(d)
    end

    function updateSize(src, event)
        s.Position = [0 0 d.Position(3:4)];
    end
end

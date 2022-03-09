function v = trimToEOSRange(v, vMinEOS, vMaxEOS, name, verbose)
%Trim value to given validity range for EOS

    % Get minimum and maximum
    vVal = value(v);
    vMin = min(vVal); 
    vMax = max(vVal);
    % Check if we are within EOS range
    below = vVal < vMinEOS;
    above = vVal > vMaxEOS;
    % Early return if we are
    if ~any(below | above), return; end
    % Issue warning if verbose output is requested
    if nargin > 4 && verbose
        if log10(vMax) > 3, fmt = 'e'; else, fmt = 'f'; end
        fmt = ['%.3', fmt, ', %.3', fmt];
        warning(['%s = [', fmt,'] out of range (', fmt,'], ', ...
            'clipping to range'], name, vMin, vMax, vMinEOS, vMaxEOS);
    end
    % Trim to range
    v(below) = vMinEOS;
    v(above) = vMaxEOS;
    
end

%{
Copyright 2009-2022 SINTEF ICT, Applied Mathematics.

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
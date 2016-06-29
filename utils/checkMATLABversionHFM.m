function varargout = checkMATLABversionHFM()
% Checks if the matlab version is greater than 2015a for compatibility with
% 'ismembertol' and 'uniquetol' functions and prints a warning if the
% matlab version is lower.
%
% SYNOPSIS:
%   v = checkMATLABversionHFM();
%   checkMATLABversionHFM();
%
% RETURNS:
%   v (OPTIONAL) - Matlab version as a character array (ex: 2015a)

%{
Copyright 2009-2015: SINTEF ICT, Applied Mathematics.

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

v = version('-release');
if str2double(v(1:4))<2015
    warning('mrst:versionIncompatibleForHFMmodule', ...
        ['The MATLAB version installed on your computer predates 2015a. \n\n'...
        'Please download MATLAB version 2015a or laster for complete '...
        'functionality with the HFM module. Lower versions of matlab '...
        'do not support the triangulation feature required to run 3D examples.']);
end

if nargout > 0, varargout{1} = v; end
end
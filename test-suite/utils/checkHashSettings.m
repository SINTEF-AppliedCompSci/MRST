function checkHashSettings()
% Utility function that checks has settings and issues warning if hashing
% is enabled.
    if ~mrstSettings('get', 'useHash')
        msg = sprintf([ ...
            'Computation of hash values is disabled in the current MRST settings.\n'         , ...
            'However, the test-suite module relies heavily on hash values, and hashing \n'   , ...
            'will be used in this example. Depending on your system, this may not work as \n', ...
            'expected. Full support for hashing will be included in a future release.'       ]);
        warning(msg); %#ok
    end
end

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
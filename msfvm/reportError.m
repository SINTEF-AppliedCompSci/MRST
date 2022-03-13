function error = reportError(reference, ms)
% Simple error helper for MsFV examples

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


    error2 = norm(reference-ms)./norm(reference);

    e = abs(reference - ms);

    M = max(e)/max(abs(reference));
    m = min(e)/min(abs(reference));

    fprintf('ERROR:\n');
    fprintf('\t2: %2.8f\n\tSup: %2.8f\n\tMinimum %2.8f\n', error2, M, m);
    error.average = error2;
    error.max = M;
    error.min = m;
end

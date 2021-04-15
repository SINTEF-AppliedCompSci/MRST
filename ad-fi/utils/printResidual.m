function printResidual(residuals, gmresits, eqnnames, iteration, CNV, MB)
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

    if iteration == 1
        fprintf('%-9s', eqnnames{:})
        fprintf('\n');
    end
    CNV = CNV(CNV ~= 0);
    MB  = MB(CNV ~= 0);
    fprintf('%8.2e ', residuals);
    fprintf('** CNV: ');
    fprintf('%2.2e ', CNV);
    fprintf('MB: ');
    fprintf('%2.2e ', MB);
    if ~isempty(gmresits)
        fprintf('** %d GMRES iterations', gmresits(2));
    end
    fprintf('\n');
end

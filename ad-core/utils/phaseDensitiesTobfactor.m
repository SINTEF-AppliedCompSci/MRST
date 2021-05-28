function b = phaseDensitiesTobfactor(rho, rhoS, dissolved)
% Convert densities to b-facctors, accounting for dissolution

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

    numPh = numel(rhoS);
    b = rho;
    
    for i = 1:numPh
        factor = rhoS(i);
        if ~isempty(dissolved)
            for j = 1:numPh
                r_ph = dissolved{j}{i};
                if ~isempty(r_ph)
                    factor = factor + rhoS(j).*r_ph;
                end
            end
        end
        if iscell(b)
            b{i} = rho{i}./factor;
        else
            b(:, i) = rho(:, i)./factor;
        end
    end
end

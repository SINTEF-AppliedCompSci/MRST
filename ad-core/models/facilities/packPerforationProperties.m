function packed = packPerforationProperties(W, p, mob, rho, dissolved, comp, wellvars, wellvarnames, varmaps, wellmap, ix)
% Extract variables corresponding to a specific well
%
% SYNOPSIS:
%  packed = packPerforationProperties(W, p, mob, rho, dissolved, comp, wellvars, wellvarnames, varmaps, wellmap, ix)
%
% REQUIRED PARAMETERS:
%  W       - Well struct for the specific well under consideration.
%
%  p       - Reservoir pressure for all cells in the domain.
%
%  mob     - Cell array with phase mobilities in each cell in the
%            reservoir. Number of active phases long.
%
%  rho     - Cell array with phase densities in each cell in the
%            reservoir. Number of active phases long.
%
%  dissolved - Black-oil specific array of rs/rv. See the function
%            ThreePhaseBlackOilModel>getDissolutionMatrix
%            for details. Should be empty for models without dissolution
%            ratios.
%
%  comp    - Cell array of components. Each entry should contain all
%            reservoir cell values of that component, subject to whatever
%            ordering is natural for the model itself.
%
%  wellvars - Extended variable set added by the FacilityModel (see
%             SimpleWell>getExtraPrimaryVariables)
%
%  wellvarnames, varmaps, wellmap, ix - Internal book-keeping for
%  additional added primary variables.
%
% RETURNS:
%   packed - Struct where the values for the perforations of the well has
%            been extracted.
%
% NOTE: This function is intended for internal use in FacilityModel.
%
% SEE ALSO:
%   FacilityModel

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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

    wc = W.cells;
    packed = struct();
    packed.pressure = p(wc);
    packed.mob = getCellSubset(mob, wc);
    packed.rho = getCellSubset(rho, wc);
    packed.dissolved = getComponentCellSubset(dissolved, wc);
    packed.components = getCellSubset(comp, wc);
    
    % Extra variables outside of standard subset
    varw = getVariableSubsetWell(wellvars, varmaps, ix);
    renum = wellmap(ix, wellmap(ix, :) > 0);
    varw = varw(renum);
    packed.extravars = varw;
    packed.extravars_names = wellvarnames;
end

function celldata = getComponentCellSubset(celldata, wc)
    for i = 1:numel(celldata)
        for j = 1:numel(celldata{i});
            if ~isempty(celldata{i}{j})
                celldata{i}{j} = celldata{i}{j}(wc);
            end
        end
    end
end

function subset = getCellSubset(celldata, wc)
    subset = cell(size(celldata));
    for i = 1:numel(subset)
        if ~isempty(celldata{i})
            subset{i} = celldata{i}(wc);
        end
    end
end

function subset = getVariableSubsetWell(vars, wellMap, ix)
    subset = cell(size(vars));
    for i = 1:numel(subset)
        subset{i} = vars{i}(wellMap{i} == ix);
    end
end
function [eqs, qBC, qRes, BCTocellMap, qSRC, srcCells, bcCells] = addFluxesFromSourcesAndBC(model, eqs, pressure, rho, mob, s, forces)
%Add in fluxes imposed by sources and face boundary conditions
%
% DESCRIPTION:
%   Utility function for updating residual conservation equations in the AD
%   framework with additional fluxes due to boundary conditions and
%   sources. Wells are handled separately in WellModel.
%
% REQUIRED PARAMETERS:
%   model      - Simulation model (subclass of ReservoirModel).
%
%   eqs        - Residual conservation of mass-equations for each phase, in
%                the order WATER, OIL, GAS (with any inactive phases
%                omitted) as a cell array.
%
%  (All the following arguments are cell arrays, with length equal to the
%  number of the active phases, with the values for each cell in each
%  entry unless otherwise noted)
%
%   pressure   - Phase pressures
%   rho        - Surface densities (one value per phase)
%   mob        - Phase mobilities
%   s          - Phase saturations
%
%   forces     - Struct containing .src and .bc fields for sources and
%                boundary conditions respectively.
%
% RETURNS:
%   eqs         - Phase conservation equations with added fluxes.
%
%   qBC         - Phase fluxes due to BC at standard conditions.
%
%   BCTocellMap - Matrix mapping qBC to cells.
%
%   qSRC        - Phase fluxes due to source terms.
%
%   srcCells    - List of cells, mapping qSRC to cells.
%
% SEE ALSO:
%   getBoundaryConditionFluxesAD, getSourceFluxesAD, addSource, addBC

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
    [qBC, qSRC, qRes] = deal(cell(numel(mob), 1));
    [BCTocellMap, srcCells, bcCells] = deal([]);
        
    if isempty(forces.bc) && isempty(forces.src)
        return
    end
    if (isprop(model, 'disgas') && model.disgas) ||...
       (isprop(model, 'vapoil') && model.vapoil)
        warning(['Boundary conditions and source terms do not fully support', ...
                 ' problems with vapoil or disgas active!']);
    end

    if ~isempty(forces.bc)
        % Setup the fluxes from the boundary condition
        %[qBC, BCTocellMap, bcCells, qRes] = getBoundaryConditionFluxesAD(model, pressure, rho, mob, s, forces.bc);
        b = cellfun(@(c1, c2) c1./c2, ...
                    rho, ...
                    mat2cell(model.getSurfaceDensities(),1, repmat(1, 1, sum(model.getActivePhases))), ...
                    'uniformoutput', false);
        [qBC, BCTocellMap, bcCells, qRes] = getBoundaryConditionFluxesAD(model, pressure, s, mob, rho, b, forces.bc);
        
        for i = 1:numel(qBC)
            % Subtract fluxes
            if isempty(eqs{i})
                continue
            end
            eqs{i}  = eqs{i} - BCTocellMap*qBC{i};
        end
    end
    
    if ~isempty(forces.src)
        % Fluxes from source terms
        [qSRC, ~, srcCells] = getSourceFluxesAD(model, mob, s, forces.src);
        for i = 1:numel(qSRC)
            if isempty(eqs{i})
                continue
            end
            % Subtract fluxes
            eqs{i}(srcCells)  = eqs{i}(srcCells) - qSRC{i};
        end
    end
end

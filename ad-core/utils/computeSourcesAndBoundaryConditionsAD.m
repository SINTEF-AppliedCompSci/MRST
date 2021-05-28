function [src, bc] = computeSourcesAndBoundaryConditionsAD(model, pressure, s, mob, rho, dissolved, forces)
%Compute phase-pseudocomponent source terms (compatible with AD codes)
%
% SYNOPSIS:
%   [src, bc] = computeSourcesAndBoundaryConditionsAD(model, pressure, s, mob, rho, dissolved, forces)
%
% REQUIRED PARAMETERS:
%   model     - Subclass of ReservoirModel for which the source terms are
%               to be computed.
%   pressure  - Reservoir pressures (cell array, one pressure per phase)
%   s         - Phase saturations (cell array, one saturation per phase)
%   mob       - Phase mobilities (cell array, one mobilit per phase)
%   rho       - Phase densities (including contributions from dissolved
%               phases. For a black-oil style model, this is
%               rhoO = bO.*(rs*rhoGS + rhoOS), and not the pseudocomponent
%               density in the phase bO.*rhoOS!
%   dissolved - Dissolution matrix for the properties. See
%               `getDissolutionMatrix` in `ThreePhaseBlackOilModel`.
%   forces    - Struct containing standard MRST driving forces.
%               Specifically, this routine uses the src and bc fields. All
%               other fields are ignored. For well source terms, see the
%               `FacilityModel`.
%
% RETURNS:
%   src,bc  - Structs for sources and boundary conditions respectively.
%             Each struct uses the same format with the following fields:
%                 'phaseMass'    - Cell array of mass source terms per
%                 phase pseudocomponent, accounting for dissolved fractions.
%                 'phaseVolume'  - Cell array of volumetric source terms
%                 per phase at reservoir conditions.
%                 'components'   - Empty cell array for inserting component
%                 source terms. Components are not the responsibility of
%                 this function, but we add the field to ensure that the
%                 structure is normalized.
%                 'mapping'      - Either empty or a matrix used to map the
%                 source terms into aggregate per-cell values. This matrix
%                 is required when multiple source terms are defined in the
%                 same block (e.g. two faces for a cell) since Matlab
%                 overwrites repeat indices instead of summing them.
%                 'sourceCells'  - List of cells the source terms should be
%                 added to.
%
% NOTE: 
%   The practical implementation of boundary conditions is normally done
%   through the gateway ReservoirModel>addBoundaryConditionsAndSources
%   routine, which uses this routine directly.
%
% SEE ALSO:
%   `equationsOilWater`, `equationsBlackOil`,
%   `ReservoirModel>addBoundaryConditionsAndSources`

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

    hasBC  = isfield(forces, 'bc')  && ~isempty(forces.bc);
    hasSRC = isfield(forces, 'src') && ~isempty(forces.src);
    rhoS = model.getSurfaceDensities();
    
    [qVolBC, qResBC, qVolSRC, b] = deal(cell(1, numel(mob)));
    [bcCells, srcCells, BCTocellMap, BCToSourceMap] = deal([]);

    if hasBC || hasSRC
        b = phaseDensitiesTobfactor(rho, rhoS, dissolved);
        if hasBC
            % Setup the fluxes from the boundary condition
            [qVolBC, BCTocellMap, bcCells, qResBC] = getBoundaryConditionFluxesAD(model, pressure, s, mob, rho, b, forces.bc);
        end

        if hasSRC
            % Fluxes from source terms
            [qVolSRC, BCToSourceMap, srcCells] = getSourceFluxesAD(model, mob, s, forces.src);
        end
    end
    src = getContributionsStruct(forces.src, qVolSRC, b, rhoS, srcCells, dissolved, BCToSourceMap);
    bc = getContributionsStruct(forces.bc, qVolBC, b, rhoS, bcCells, dissolved, BCTocellMap, qResBC);
end

function src = getContributionsStruct(force, q_s, b, rhoS, cells, dissolved, map, q_r)
    nPh = numel(q_s);
    if nargin < 8
        q_r = q_s;
        for i = 1:nPh
            q_r{i} = q_s{i}./b{i}(cells);
        end
    end
    
    if ~isempty(dissolved) && ~isempty(force)
        q_s0 = q_s;
        q_t = 0;
        for i = 1:numel(q_s)
            q_t = q_t + value(q_s{i});
        end
        isInj = q_t > 0;
        for i = 1:nPh
            for j = 1:nPh
                % Add dissolution of component i into phase j
                r_cell = dissolved{i}{j};
                if isempty(r_cell)
                    continue
                end
                
                if isfield(force, 'dissolution')
                    % Note: If dissolution is specified for rate
                    % sources/bc, the total injected mass will change as
                    % the dissolved mass is not accounted for. This feature
                    % is primarily intended for pressure bc.
                    ds_mat = force.dissolution;
                    if ndims(ds_mat) == 3
                        % One value per cell, per phase and component
                        r_inj = ds_mat(:, i, j);
                    else
                        % One value for all cells, per phase and component
                        r_inj = ds_mat(i, j);
                    end
                else
                    r_inj = 0;
                end
                r = isInj.*r_inj + ~isInj.*r_cell(cells);
                q_s{i} = q_s{i} + r.*q_s0{j};
            end
        end
    end
    if nargin > 6
        mm = map(cells, :);
%         for i = 1:numel(q_s)
%             q_s{i} = mm*q_s{i};
%             q_r{i} = mm*q_r{i};
%         end
    else
        mm = [];
    end
    
    srcMass = q_s;
    for i = 1:nPh
        srcMass{i} = srcMass{i}.*rhoS(i);
    end
    src = struct('phaseMass',   {srcMass}, ...
                 'phaseVolume', {q_r}, ...
                 'components',  {[]}, ...
                 'mapping',     mm, ...
                 'sourceCells', cells);
end

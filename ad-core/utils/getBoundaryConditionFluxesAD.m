function [qSurf, BCTocellMap, BCcells, qRes] = getBoundaryConditionFluxesAD(model, pressure, s, mob, rho, b, bc)
%Get boundary condition fluxes for a given set of values
%
% SYNOPSIS:
%   [qSurf, BCTocellMap, BCcells] = getBoundaryConditionFluxesAD(model, pressure, s, mob, rho, b, bc)
%
% DESCRIPTION:
%   Given a set of boundary conditions, this function computes the fluxes
%   induced for a given set of reservoir parameters (density, mobility,
%   saturations etc).
%
% REQUIRED PARAMETERS:
%
%   model      - Subclass of ReservoirModel implementing the current
%                simulation model.
%
%   pressure   - Cell values of pressure. Should be a nph long cell array,
%                containing the phase pressures.
%
%   rho        - Surface densities of each phase, as a nph long cell array.
%
%   s          - Phase saturations per cell, as a nph long array.
%
%   bc         - Boundary condition struct, with valid .sat field with
%                length nph. Typically made using addBC, pside or fluxside.
%
% RETURNS:
%   qSurf      - Cell array of phase fluxes.
%
%   BCTocellMap - Matrix used to add in bc fluxes to cells. Implemented as
%                 a matrix to efficiently account for cells with multiple
%                 faces with boundary conditions.
%
%   cells       - The cells affected by boundary conditions.
%
% SEE ALSO:
%   addBC, pside, fluxside

%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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

% Basic quanitites
T = model.operators.T_all(bc.face);
G = model.G;
nPh = sum(model.getActivePhases);
N = G.faces.neighbors(bc.face,:);

% Validation
assert(size(bc.sat, 2) == nPh, ...
    ['Wrong number of columns in BC sat field: Expected columns', ...
    num2str(nPh), ', but input had ', num2str(size(bc.sat, 2)), ' columns.']);
assert(~any(all(N > 0, 2)),'bc on internal boundary');

hasOutsideMob = isfield(bc, 'mob');
hasOutsideRho = isfield(bc, 'rho');

if ~hasOutsideMob
    ss = sum(bc.sat, 2);
    % Values should either sum to zero or one
    assert(all(ss - 1 < sqrt(eps) | ss < sqrt(eps)));
end
% Mapping
BCcells = sum(N, 2);
nbc = numel(bc.face);
% This mapping takes us from cell values to bc values and vice versa. We
% use sparse matrices to add in fluxes because using sub-indices will
% overwrite values when multiple BC are applied to the same cell (i.e. a
% corner cell with BC on multiple sides).
cellToBCMap = sparse((1:nbc)', BCcells, 1, nbc, G.cells.num);
BCTocellMap = cellToBCMap';

% Gravity gradient per bc face 
if any(strcmpi(G.type, 'topSurfaceGrid'))
   dzbc = model.gravity(3) * (G.cells.z(BCcells) - G.faces.z(bc.face));
else
   g = model.getGravityVector();
   dz = G.cells.centroids(BCcells, :) - G.faces.centroids(bc.face,:);
   dzbc = -dz*g';
end

isP = reshape(strcmpi(bc.type, 'pressure'), [], 1);
isSF = reshape(strcmpi(bc.type, 'flux'), [], 1);
isRF = ~(isP | isSF);

[qSurf, qRes] = deal(cell(nPh,1));

% Use sat field to determine what any inflow cells produce.
sat = bc.sat;
noSat = all(sat == 0, 2);
hasNoSat = any(noSat);

% Store total mobility
totMob = zeros(size(sat, 1), 1);
for i = 1:nPh
    totMob = totMob + cellToBCMap*mob{i};
end

rhoS = model.getSurfaceDensities();
for i = 1:nPh
    if isa(totMob, 'ADI')
        sample = totMob;
    else
        sample = pressure{1};
    end
    zeroAD = model.AutoDiffBackend.convertToAD(zeros(nbc, 1), sample);
    [q_s, q_r] = deal(zeroAD);
    
    pBC   = cellToBCMap*pressure{i};
    bBC_in = cellToBCMap*b{i};
    rhoBC_in = cellToBCMap*rho{i};
    mobBC_in = cellToBCMap*mob{i};
    
    if hasOutsideRho
        rhoBC_out = bc.rho(:, i);
        bBC_out = rhoBC_out./rhoS(i);
    else
        rhoBC_out = rhoBC_in;
        bBC_out = bBC_in;
    end
    if hasOutsideMob
        mobBC_out = bc.mob(:, i);
    else
        mobBC_out = mobBC_in;
    end
    
    sBC   = cellToBCMap*s{i};
    
    if hasNoSat
        % If no saturations are defined, we explicitly set it to mirror the
        % cell values on the other side of the interface
        sBC = double(sBC);
        sat(noSat, i) = sBC(noSat);
    end
    
    if any(isP)
        % Treat pressure BC
        if isa(model, 'ThreePhaseCompositionalModel') && i > model.water
            sT = (sBC(isP) + sat(isP, i));
            sT(double(sT) == 0) = 1e-8;
            rhoF = (sBC(isP).*rhoBC_in(isP) + sat(isP, i).*rhoBC_out(isP))./sT;
        else
            rhoF = (rhoBC_in(isP) + rhoBC_out(isP))./2;
        end
        dP = bc.value(isP) - pBC(isP) - rhoF.*dzbc(isP);
        
        % Determine if pressure bc are injecting or producing
        injDir = dP > 0;

        injP = isP;
        injP(isP) = injDir;

        if any(~injDir)
            % Write out the flux equation over the interface
            subs = isP & ~injP;
            q_res = mobBC_in(subs).*T(subs).*dP(~injDir);
            q_s(subs) = bBC_in(subs).*q_res;
            q_r(subs) = q_res;
            clear subs
        end

        if any(injDir)
            % In this case, pressure drives flow inwards, we get the injection rate
            % determined by the sat field
            subs = isP & injP;
            if hasOutsideMob
                q_res = mobBC_out(subs).*T(subs).*dP(injDir);
            else
                q_res = totMob(subs).*T(subs).*dP(injDir).*sat(subs, i);
            end
            q_s(subs) = bBC_out(subs).*q_res;
            q_r(subs) = q_res;
            clear subs
        end
    end
    % Treat flux / Neumann BC
    
    % ------ Fluxes given at surface conditions ----- %
    injNeu = bc.value > 0;
    
    subs = isSF &  injNeu;
    % Injection
    if any(subs)
        q_s(subs) = bc.value(subs).*sat(subs, i);
        q_r(subs) = bc.value(subs).*sat(subs, i)./bBC_in(subs);
    end
    
    subs = isSF & ~injNeu;
    % Production
    if any(subs)
        % Production fluxes, use fractional flow of total mobility to
        % estimate how much mass will be removed.
        f = mobBC_in(subs)./totMob(subs);
        tmp = f.*bc.value(subs);
        q_s(subs) = tmp;
        q_r(subs) = tmp./bBC_in(subs);
    end
    % ------ Fluxes given at reservoir conditions ----- %
    subs = isRF &  injNeu;
    % Injection
    if any(subs)
        tmp = bc.value(subs).*sat(subs, i);
        q_s(subs) = tmp.*bBC_in(subs);
        q_r(subs) = tmp;
    end
    subs = isRF & ~injNeu;
    % Production
    if any(subs)
        f = mobBC_in(subs)./totMob(subs);
        tmp = f.*bc.value(subs);
        q_s(subs) = tmp.*bBC_in(subs);
        q_r(subs) = tmp;
    end

    qSurf{i} = q_s;
    qRes{i} = q_r;
end
end


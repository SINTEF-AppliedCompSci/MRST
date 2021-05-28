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
%   b          - Inverse formation volume factors.
%
%   bc         - Boundary condition struct, with valid .sat field with
%                length nph. Typically made using addBC, pside or fluxside.
%
% RETURNS:
%   qSurf       - Cell array of phase fluxes at surface conditions.
%
%   BCTocellMap - Matrix used to add in bc fluxes to cells. Implemented as
%                 a matrix to efficiently account for cells with multiple
%                 faces with boundary conditions.
%
%   BCcells     - The cells affected by boundary conditions.
%
%   qRes        - Cell array of phase fluxes at reservoir conditions.
%
% SEE ALSO:
%   addBC, pside, fluxside

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

% Basic quanitites
T = model.operators.T_all(bc.face);
G = model.G;
nph = sum(model.getActivePhases);
N = G.faces.neighbors(bc.face,:);

% Validation
assert(size(bc.sat, 2) == nph, ...
    ['Wrong number of columns in BC sat field: Expected ', num2str(nph), ...
     ' columns, but input had ', num2str(size(bc.sat, 2)), ' columns.']);
assert(~any(all(N > 0, 2)),'bc on internal boundary');

hasOutsideMob = isfield(bc, 'mob');
hasOutsideRho = isfield(bc, 'rho');

if ~hasOutsideMob
    ss = sum(bc.sat, 2);
    % Values should either sum to zero or one
    assert(all(ss - 1 < sqrt(eps) | ss < sqrt(eps)), ...
        'Boundary conditions should have ''sat'' field which has row-sum 1.');
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

[qSurf, qRes] = deal(cell(nph,1));

% Use sat field to determine what any inflow cells produce.
sat = bc.sat;
noSat = all(sat == 0, 2);
hasNoSat = any(noSat);

% Store total mobility

isCompositional = isa(model, 'ThreePhaseCompositionalModel');
isTransport = isa(model, 'TransportNaturalVariablesModel');

rhoS = model.getSurfaceDensities();

[pressureF, rhoF, bF, mobF, sF] = deal(cell(nph, 2));
totMob = zeros(nbc, 1);
for i = 1:nph
    % First column is outside, second is inside

    % Store pressure on inside and outside
    p_inside   = cellToBCMap*pressure{i};
    p_bnd = p_inside;
    if size(bc.value, 2) > 1
        p_bnd(isP, i) = bc.value(isP, i);
    else
        p_bnd(isP) = bc.value(isP);
    end
    pressureF{i, 1} = p_bnd;
    pressureF{i, 2} = p_inside;
    % Density and b-factors
    bBC_inside = cellToBCMap*b{i};
    rhoBC_inside = cellToBCMap*rho{i};
    if hasOutsideRho
        rhoBC_bnd = bc.rho(:, i);
        bBC_bnd = rhoBC_bnd./rhoS(i);
    else
        rhoBC_bnd = rhoBC_inside;
        bBC_bnd = bBC_inside;
    end
    rhoF{i, 1} = rhoBC_bnd;
    rhoF{i, 2} = rhoBC_inside;
    bF{i, 1} = bBC_bnd;
    bF{i, 2} = bBC_inside;
    % Mobility
    mob_inside = cellToBCMap*mob{i};
    if hasOutsideMob
        mob_bnd = bc.mob(:, i);
    else
        mob_bnd = mob_inside;
    end
    mobF{i, 1} = mob_bnd;
    mobF{i, 2} = mob_inside;
    totMob = totMob + mob_inside;
    % Saturation
    s_inside   = cellToBCMap*s{i};
    if hasNoSat
        % If no saturations are defined, we explicitly set it to mirror the
        % cell values on the other side of the interface
        s_inside = value(s_inside);
        sat(noSat, i) = s_inside(noSat);
    end
    sF{i, 1} = sat(:, i);
    sF{i, 2} = s_inside;
end

if ~hasOutsideMob
    for i = 1:nph
        mobF{i, 1} = totMob.*sF{i, 1};
    end
end

if isTransport && isCompositional
    sT = cell(1, 2);
    sT{2} = zeros(numel(T), 1);
    for i = 1:nph
        sT{2} = sT{2} + sF{i, 2};
    end
    sT{1} = sum(sat, 2);
end

rhoAvgF = cell(nph, 1);
for i = 1:nph
    if isCompositional && i > model.water
        if isTransport
            sL = sF{i, 1}./sT{1};
            sR = sF{i, 2}./sT{2};
        else
            sL = sF{i, 1};
            sR = sF{i, 2};
        end
        sTf = sL + sR;
        sTf(value(sTf) == 0) = 1e-8;
        rhoAvgF{i} = (sL.*rhoF{i, 1} + sR.*rhoF{i, 2})./sTf;
    else
        rhoAvgF{i} = (rhoF{i, 1} + rhoF{i, 2})./2;
    end
end

if any(isRF)
    mrstModule add sequential
    G = cell(1, nph);

    mobC = cell(1, nph);
    for i = 1:nph
        G{i} = dzbc.*rhoAvgF{i};
        mobC{i} = vertcat(mobF{i, 2}, mobF{i, 1});
    end
    vT = sum(bc.value, 2);

    nf = numel(vT);
    if 0
        upstr = @(flag, v) flag.*v(1:nf, :) + ~flag.*v(nf+1:end, :);
        q_ph = computeSequentialFluxes([], G, vT, T, mobC, {}, {}, upstr, 'potential');
    else
        upstr = @(flag, v) flag.*v(1:nf, :) + ~flag.*v(nf+1:end, :);
        q_ph = computeSequentialFluxes([], G, -vT, T, mobC, {}, {}, upstr, 'potential');
        for i = 1:numel(q_ph)
            q_ph{i} = - q_ph{i};
        end
    end
end

for i = 1:nph
    if size(bc.value, 2) == 1
        bc_v = bc.value;
    else
        assert(size(bc.value, 2) == nph, ...
        'Boundary conditions should have either one value or one value per phase, for each face');
        bc_v = bc.value(:, i);
    end

    if isa(totMob, 'ADI')
        sample = totMob;
    else
        sample = pressure{1};
    end
    zeroAD = model.AutoDiffBackend.convertToAD(zeros(nbc, 1), sample);
    [q_s, q_r] = deal(zeroAD);

    if any(isP)
        % Treat pressure BC
        rhoF = rhoAvgF{i}(isP);
        dP = pressureF{i, 1}(isP) - pressureF{i, 2}(isP) - rhoF.*dzbc(isP);

        % Determine if pressure bc are injecting or producing
        injDir = dP > 0;

        injP = isP;
        injP(isP) = injDir;

        if any(~injDir)
            % Write out the flux equation over the interface
            subs = isP & ~injP;
            q_res = mobF{i, 2}(subs).*T(subs).*dP(~injDir);
            q_s(subs) = bF{i, 2}(subs).*q_res;
            q_r(subs) = q_res;
            clear subs
        end

        if any(injDir)
            % In this case, pressure drives flow inwards, we get the injection rate
            % determined by the sat field
            subs = isP & injP;
            if hasOutsideMob
                q_res = mobF{i, 1}(subs).*T(subs).*dP(injDir);
            else
                q_res = totMob(subs).*T(subs).*dP(injDir).*sat(subs, i);
            end
            q_s(subs) = bF{i, 1}(subs).*q_res;
            q_r(subs) = q_res;
            clear subs
        end
    end
    % Treat flux / Neumann BC

    % ------ Fluxes given at surface conditions ----- %
    injNeu = bc_v > 0;
    if any(isSF)
        subs = isSF &  injNeu;
        % Injection
        if any(subs)
            q_s(subs) = bc_v(subs).*sat(subs, i);
            q_r(subs) = bc_v(subs).*sat(subs, i)./bF{i, 1}(subs);
        end

        subs = isSF & ~injNeu;
        % Production
        if any(subs)
            % Production fluxes, use fractional flow of total mobility to
            % estimate how much mass will be removed.
            f = mobF{i, 2}(subs)./totMob(subs);
            tmp = f.*bc_v(subs);
            q_s(subs) = tmp;
            q_r(subs) = tmp./bF{i, 2}(subs);
        end
    end
    % ------ Fluxes given at reservoir conditions ----- %
    if any(isRF)
        injNeuR = q_ph{i} > 0;
        subs = isRF &  injNeuR;
        % Injection
        if any(subs)
            tmp = q_ph{i}(subs);
            q_s(subs) = tmp.*bF{i, 1}(subs);
            q_r(subs) = tmp;
        end
        subs = isRF & ~injNeuR;
        % Production
        if any(subs)
            tmp = q_ph{i}(subs);
            q_s(subs) = tmp.*bF{i, 2}(subs);
            q_r(subs) = tmp;
        end
    end
    qSurf{i} = q_s;
    qRes{i} = q_r;
end
end

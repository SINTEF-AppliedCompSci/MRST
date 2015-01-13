function [qRes, BCTocellMap, BCcells] = getBoundaryConditionFluxesAD(model, pressure, rho, mob, b, s, bc)

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
ss = sum(bc.sat, 2);
assert(all(ss == 1 | ss == 0));

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
g = model.getGravityVector();
dz = G.cells.centroids(BCcells, :) - G.faces.centroids(bc.face,:);
dzbc = dz*g';

isP = reshape(strcmpi(bc.type, 'pressure'), [], 1);

qRes = cell(nPh,1);

% Use sat field to determine what any inflow cells produce.
sat = bc.sat;
noSat = all(sat == 0, 2);
hasNoSat = any(noSat);

% Store total mobility
totMob = double2ADI(zeros(nbc, 1), mob{1});
for i = 1:nPh
    totMob = totMob + cellToBCMap*mob{i};
end

for i = 1:nPh
    
    q = double2ADI(zeros(nbc, 1), mob{i});
    
    pBC   = cellToBCMap*pressure{i};
    rhoBC = cellToBCMap*rho{i};
    bBC   = cellToBCMap*b{i};
    mobBC = cellToBCMap*mob{i};
    sBC   = cellToBCMap*s{i};
    
    if hasNoSat
        % If no saturations are defined, we explicitly set it to mirror the
        % cell values on the other side of the interface
        sBC = double(sBC);
        sat(noSat, i) = sBC(noSat);
    end
    
    % Treat pressure BC
    dP = bc.value(isP) - pBC(isP) + rhoBC(isP).*dzbc(isP);

    % Determine if pressure bc are injecting or producing
    injDir = dP > 0;
    
    injP = isP;
    injP(isP) = injDir;
    
    if any(~injDir)
        % Write out the flux equation over the interface
        q(isP & ~injP)  = bBC(isP & ~injP).*mobBC(isP & ~injP).*T(isP & ~injP).*dP(~injDir);
    end
    
    if any(injDir)
        % In this case, pressure drives flow inwards, we get the injection rate
        % determined by the sat field
        q( isP & injP)  = bBC( isP & injP).*totMob(isP & injP).*T( isP & injP).*dP(injDir).*sat(isP & injP, i);
    end
    % Treat flux / Neumann BC
    injNeu = bc.value > 0;
    if any(~isP &  injNeu)
        % Injection
        q(~isP &  injNeu) = bc.value(~isP & injNeu).*sat(~isP & injNeu, i);
    end
    if any(~isP & ~injNeu)
        % Production
        q(~isP & ~injNeu) = bc.value(~isP & ~injNeu).*sBC(~isP & ~injNeu);
    end
    qRes{i} = q;
end
end
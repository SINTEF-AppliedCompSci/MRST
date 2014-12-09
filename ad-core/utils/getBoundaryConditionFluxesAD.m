function [qRes, bc2cellMap, bc_cell] = getBoundaryConditionFluxesAD(model, gdz, pressure, rho, mob, b, s, bc)

% Basic quanitites
T = model.operators.T_all(bc.face);
G = model.G;
nPh = sum(model.getActivePhases);
N = G.faces.neighbors(bc.face,:);

% Validation
assert(size(bc.sat, 2) == nPh);
assert(~any(all(N > 0, 2)),'bc on internal boundary');

% Mapping
bc_cell = sum(N, 2);
nbc = numel(bc.face);
cell2bcMap = sparse((1:nbc)', bc_cell, 1, nbc, G.cells.num);
bc2cellMap = cell2bcMap';

dzbc = gdz(bc.face);



isP = reshape(strcmpi(bc.type, 'pressure'), [], 1);

qRes = cell(numel(pressure),1);

% Use sat field to determine what any inflow cells produce.
sat = bc.sat;
noSat = all(sat == 0, 2);
hasNoSat = any(noSat);

% Store total mobility
totMob = double2ADI(zeros(nbc, 1), mob{1});
for i = 1:nPh
    totMob = totMob + cell2bcMap*mob{i};
end

for i = 1:nPh
    
    q = double2ADI(zeros(nbc, 1), mob{i});
    
    pBC   = cell2bcMap*pressure{i};
    rhoBC = cell2bcMap*rho{i};
    bBC   = cell2bcMap*b{i};
    mobBC = cell2bcMap*mob{i};
    sBC   = cell2bcMap*s{i};
    
    if hasNoSat
        % If no saturations are defined, we explicitly set it to mirror the
        % cell values on the other side of the interface
        sBC = double(sBC);
        sat(noSat, i) = sBC(noSat);
    end
    
    % Treat pressure BC
    dP = bc.value(isP) - pBC(isP) + rhoBC(isP).*dzbc(isP);

    % Determine if pressure bc are injecting or producing
    inj = dP > 0;
    
    injP = isP;
    injP(isP) = inj;
    
    % Write out the flux equation over the interface
    q(isP & ~injP)  = -bBC(isP & ~injP).*mobBC(isP & ~injP).*T(isP & ~injP).*dP(~inj);

    % In this case, pressure drives flow inwards, we get the injection rate
    % determined by the sat field
    q( isP & injP)  = -bBC( isP & injP).*totMob(isP & injP).*T( isP & injP).*dP(inj).*sat(isP & injP, i);

    % Treat flux / Neumann BC
    inj = bc.value > 0;
    % Injection
    q(~isP &  inj) = -bc.value(~isP & inj).*sat(~isP & inj, i);
    % Production
    q(~isP & ~inj) = -bc.value(~isP & ~inj).*sBC(~isP & ~inj);
    
    qRes{i} = q;
end
end
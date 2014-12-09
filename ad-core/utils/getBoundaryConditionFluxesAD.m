function [qRes, bc2cellMap, bc_cell] = getBoundaryConditionFluxesAD(model, gdz, pressure, rho, mob, b, s, bc)

% Basic quanitites
T = model.operators.T_all(bc.face);
G = model.G;
nPh = sum(model.getActivePhases);
N = G.faces.neighbors(bc.face,:);

% Validation
assert(all(strcmp(bc.type,'pressure')), 'only pressure bc allowed');
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
    
    % Treat pressure BC
    dP = bc.value(isP) - pBC(isP) + rhoBC(isP).*dzbc(isP);

    % Determine if pressure bc are injecting or producing
    inj = dP > 0;
    injP = inj & isP;

    q(~injP)  = -bBC(~injP).*mobBC(~injP).*T(~injP).*dP(~inj);

    % Pressure drives flow inwards, we get the injection rate
    % determined by the sat field
    q( injP)  = -bBC( injP).*totMob(injP).*T( injP).*dP(inj).*bc.sat(injP, i);

    % Treat flux / Neumann BC
    q(~isP) = bc.value(~isP).*bc.sat(~isP, i);

    qRes{i} = q;
end
end
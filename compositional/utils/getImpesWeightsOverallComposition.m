function w = getImpesWeightsOverallComposition(model, state, dt)
state = model.EOSModel.updateAfterConvergence(state, state, dt, struct());

[ncell, ncomp] = size(state.components);
s = model.operators;
[p, z, temp] = model.getProps(state, ...
    'pressure', 'components', 'T');
z = expandMatrixToCell(z);

[p, z{1:end-1}] = initVariablesADI(...
 p, z{1:end-1});

z{end} = 1;
for i = 1:ncomp-1
    z{end} = z{end} - z{i};
end
acc = cell(1, ncomp);
[xM, yM, sO, sG, rhoO, rhoG] = model.computeTwoPhaseFlowProps(state, p, temp, z);
for i = 1:ncomp
    ci = (s.pv/dt).*(rhoO.*sO.*xM{i} + rhoG.*sG.*yM{i} );
    acc{i} = ci;
end

ok = ~cellfun(@isempty, acc);
e = vertcat(acc{ok});
c = cat(e);
A = c.jac{1};
ndof = ncell*ncomp;

b = zeros(ndof, 1);
b(1:ncell) = 1/barsa;

Ap = A';
w = Ap\b;
w = reshape(w, [], ncomp);
end

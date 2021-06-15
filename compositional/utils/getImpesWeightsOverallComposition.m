function w = getImpesWeightsOverallComposition(model, state, dt)
%Undocumented Utility Function

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

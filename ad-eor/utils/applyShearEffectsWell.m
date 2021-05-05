function mob = applyShearEffectsWell(mob, q_ph, prop, model, state)
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

wIx = 1;
q_w = value(q_ph{wIx});

[b, cp] = model.getProps(state, 'ShrinkageFactors', 'polymer');

map = prop.getEvaluatedDependencies(state, 'FacilityWellMapping');
W = map.W;
wc = map.cells;

fluid = model.fluid;
s    = model.operators;
G = model.G;

bwW = value(b{wIx}(wc));
muWeffMult = model.getProp(state, 'PolymerEffViscMult');
muWMultW = value(muWeffMult(wc));


[~, wciPoly, iInxW] = getWellPolymer(W);
cpw = value(cp(wc));
% for injectors, we assume the polymer is full mixed
muWMultW(iInxW) = fluid.muWMult(cpw(iInxW));

% Maybe should also apply this for PRODUCTION wells.
muWMultW((iInxW(wciPoly==0))) = 1;

poro =  s.pv./G.cells.volumes;
poroW = poro(wc);


% the thickness of the well perforations in the cell
% TODO: the following can be put into a function
welldir = { W.dir };
i = cellfun('prodofsize', welldir) == 1;
welldir(i) = arrayfun(@(w) repmat(w.dir, [ numel(w.cells), 1 ]), ...
                      W(i), 'UniformOutput', false);
welldir = vertcat(welldir{:});
[dx, dy, dz] = cellDims(G, wc);
thicknessWell = dz;
thicknessWell(welldir == 'Y') = dy(welldir == 'Y');
thicknessWell(welldir == 'X') = dx(welldir == 'X');

if ~isfield(W, 'rR')
     error('The representative radius of the well is not initialized');
end
rR = vertcat(W.rR);


% TODO: probably it will be easier just to calculate the lateral area of
% the well bore to save some computation
% this is where the BUG of E100 is. 
VwW = bwW .* bwW.*q_w./(poroW .* rR .* thicknessWell * 2 * pi);
VwW(value(q_w) == 0) = 0;
% VwW = q_w./(poroW .* rR .* thicknessWell * 2 * pi);
shearMultW = computeShearMult(fluid, abs(VwW), muWMultW);

mob{wIx} = mob{wIx} ./ shearMultW;

end


% TODO: this function is a copy, possibly we should make it public
function [wPoly, wciPoly, iInxW] = getWellPolymer(W)
    if isempty(W)
        wPoly = [];
        wciPoly = [];
        iInxW = [];
        return
    end
    inj   = vertcat(W.sign)==1;
    polInj = cellfun(@(x)~isempty(x), {W(inj).cp});
    wPoly = zeros(nnz(inj), 1);
    W_inj = W(inj);
    wPoly(polInj) = vertcat(W_inj(polInj).cp);
    wciPoly = rldecode(wPoly, cellfun(@numel, {W_inj.cells}));

    % Injection cells
    nPerf = cellfun(@numel, {W.cells})';
    nw    = numel(W);
    perf2well = rldecode((1:nw)', nPerf);
    compi = vertcat(W.compi);
    iInx  = rldecode(inj, nPerf);
    iInx  = find(iInx);
    iInxW = iInx(compi(perf2well(iInx),1)==1);
end

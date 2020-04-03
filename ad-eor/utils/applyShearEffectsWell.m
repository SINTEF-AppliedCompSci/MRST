function shearMultW = applyShearEffectsWell(q_ph, prop, model, state)

% TODO: to be improved to avoid wasting computation
check = @(prop) isprop(model, prop) && model.(prop);

if ~(check('usingShear') || check('usingShearLog') || check('usingShearLogshrate'))
    shearMultW = 1.;
    return;
end


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
% VwW = q_w./(poroW .* rR .* thicknessWell * 2 * pi);
shearMultW = computeShearMult(fluid, abs(VwW), muWMultW);

end


% TODO: this function is a copy
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
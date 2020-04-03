function v = applyShearEffects(v, model, state)

check = @(prop) isprop(model, prop) && model.(prop);

if ~(check('usingShear') || check('usingShearLog') || check('usingShearLogshrate'))
    return;
end

% TODO: I believe there is some repeated calculation in the following
s    = model.operators;
G = model.G;

% TODO: only handle model.usingShear first
% the indice of water phase, and also the water component
wIx = 1;
vW = v{wIx, wIx};
% we discard the derivative here to save computation
vW = value(vW);

upcw = model.getProp(state, 'PhaseUpwindFlag');
upcw = upcw{wIx};

rho = model.getProp(state, 'Density');
rhow = rho{wIx}.val;
rhowf = s.faceUpstr(upcw, rhow);


poro =  s.pv./G.cells.volumes;
poroFace = s.faceAvg(poro);
faceA = G.faces.areas(s.internalConn);

% this is the water velocity based on the flux rate
Vw = vW./rhowf./(poroFace .* faceA);

muWeffMult = value(model.getProp(state, 'PolymerEffViscMult'));
muWMultf = s.faceUpstr(upcw, muWeffMult);

% TODO: we should try to obtain the derivative, potentially for better
% convergence
shearMultf = computeShearMult(model.fluid, abs(Vw), muWMultf);

ncomp = model.getNumberOfComponents;
for c = 1:ncomp
    if ~isempty(v{c, wIx})
        v{c, wIx} = v{c, wIx} ./ shearMultf;
    end
end

end
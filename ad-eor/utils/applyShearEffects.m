function v = applyShearEffects(v, model, state)
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
rhow = value(rho{wIx});
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

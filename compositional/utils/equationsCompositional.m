function [problem, state] = equationsCompositional(state0, state, model, dt, drivingForces, varargin)
opt = struct('Verbose',     mrstVerbose,...
            'reverseMode', false,...
            'resOnly',     false,...
            'iteration',   -1);

opt = merge_options(opt, varargin{:});

% Shorter names for some commonly used parts of the model and forces.
s = model.operators;
f = model.fluid;
W = drivingForces.W;

fluid = model.fluid;
compFluid = model.EOSModel.fluid;

state = model.computeFlash(state, dt, opt.iteration);
% Properties at current timestep
[p, sW, z, temp, wellSol] = model.getProps(state, ...
    'pressure', 'water', 'z', 'T', 'wellSol');
assert(all(p>0), 'Pressure must be positive for compositional model');

[p0, sW0, z0, temp0] = model.getProps(state0, ...
    'pressure', 'water', 'z', 'T');

bhp    = vertcat(wellSol.bhp);
qWs    = vertcat(wellSol.qWs);
qOs    = vertcat(wellSol.qOs);
qGs    = vertcat(wellSol.qGs);

ncomp = numel(z);
cnames = model.EOSModel.fluid.names;
if model.water
    [p, z{1:ncomp-1}, sW, qWs, qOs, qGs, bhp] = initVariablesADI(p, z{1:ncomp-1}, sW, qWs, qOs, qGs, bhp);
    primaryVars = {'pressure', cnames{1:end-1}, 'sW', 'qWs', 'qOs', 'qGs', 'bhp'};
else
    [p, z{1:ncomp-1}, qOs, qGs, bhp] = initVariablesADI(p, z{1:ncomp-1}, qOs, qGs, bhp);
    primaryVars = {'pressure', cnames{1:end-1}, 'qOs', 'qGs', 'bhp'};
end
z{end} = 1;
for i = 1:(ncomp-1)
    z{end} = z{end} - z{i};
end

[xM,  yM,  sO,  sG,  rhoO,  rhoG, muO, muG] = model.computeTwoPhaseFlowProps(state, p, temp, z);
[xM0, yM0, sO0, sG0, rhoO0, rhoG0] = model.computeTwoPhaseFlowProps(state0, p0, temp0, z0);

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

if model.water
    sO = (1-sW).*sO;
    sG = (1-sW).*sG;
    sO0 = (1-sW0).*sO0;
    sG0 = (1-sW0).*sG0;
    
    [krW, krO, krG] = model.evaluateRelPerm({sW, sO, sG});
    krW = mobMult.*krW;
else
    [krO, krG] = model.evaluateRelPerm({sO, sG});
end

krO = mobMult.*krO;
krG = mobMult.*krG;

% Compute transmissibility
T = s.T.*transMult;

% Gravity gradient per face
gdz = model.getGravityGradient();

% Oil flux
rhoOf  = s.faceAvg(rhoO);
mobO   = krO./muO;
dpO    = s.Grad(p) - rhoOf.*gdz;
upco  = (double(dpO)<=0);
vO = -s.faceUpstr(upco, mobO).*T.*dpO;

% Gas flux
rhoGf  = s.faceAvg(rhoG);
mobG   = krG./muG;
dpG    = s.Grad(p) - rhoGf.*gdz;
upcg  = (double(dpG)<=0);
vG = -s.faceUpstr(upcg, mobG).*T.*dpG;

rOvO = s.faceUpstr(upco, rhoO).*vO;
rGvG = s.faceUpstr(upcg, rhoG).*vG;

bO = rhoO./fluid.rhoOS;
bG = rhoG./fluid.rhoGS;
% EQUATIONS -----------------------------------------------------------
if model.water
    % Water flux
    muW = f.muW(p);
    bW     = fluid.bW(p);
    rhoW   = bW.*fluid.rhoWS;
    rhoW0 = fluid.bW(p0).*fluid.rhoWS;

    rhoWf  = s.faceAvg(rhoW);
    mobW   = krW./muW;
    dpW    = s.Grad(p) - rhoWf.*gdz;
    upcw  = (double(dpW)<=0);
    vW = -s.faceUpstr(upcw, mobW).*T.*dpW;
    rWvW = s.faceUpstr(upcw, rhoW).*vW;
    water = (s.pv/dt).*( rhoW.*pvMult.*sW - rhoW0.*pvMult0.*sW0 ) + s.Div(rWvW);
else
    [rWvW, upcw, bW] = deal([]);
end
if model.outputFluxes
    state = model.storeFluxes(state, rWvW, rOvO, rGvG);
end
if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, bG);
    state = model.storeMobilities(state, mobW, mobO, []);
    state = model.storeUpstreamIndices(state, upcw, upco, upcg);
end

% water equation + n component equations
[eqs, types, names] = deal(cell(1, ncomp + model.water));

if model.water
    eqs{1} = water;
    names{1} = 'water';
    types{1} = 'cell';
    woffset = 1;
else
    woffset = 0;
end

for i = 1:ncomp
    names{i+woffset} = compFluid.names{i};
    types{i+woffset} = 'cell';

    eqs{i+woffset} = (s.pv/dt).*( ...
                    rhoO.*pvMult.*sO.*xM{i} - rhoO0.*pvMult0.*sO0.*xM0{i} + ...
                    rhoG.*pvMult.*sG.*yM{i} - rhoG0.*pvMult0.*sG0.*yM0{i}) ...
          + s.Div(rOvO.*s.faceUpstr(upco, xM{i}) + rGvG.*s.faceUpstr(upcg, yM{i}));
end


% Finally, add in and setup well equations
if ~isempty(W)
    wm = model.wellmodel;
    
    offset = ncomp + woffset;
    if ~opt.reverseMode
        % Store cell wise well variables in cell arrays and send to ewll
        % model to get the fluxes and well control equations.
        wc    = vertcat(W.cells);
        pw    = p(wc);
        w_comp = vertcat(W.components);
        perf2well = getPerforationToWellMapping(W);
        a = w_comp(perf2well, :).*repmat(compFluid.molarMass, numel(wc), 1);
        w_comp = bsxfun(@rdivide, a, sum(a, 2));
        
        x_comp = cellfun(@(v) v(wc), xM, 'UniformOutput', false);
        y_comp = cellfun(@(v) v(wc), yM, 'UniformOutput', false);
        
        if model.water
            rhows = [f.rhoWS, f.rhoOS, f.rhoGS];
            bw    = {bW(wc), bO(wc), bG(wc)};
            mw    = {mobW(wc), mobO(wc), mobG(wc)};
            sat   = {sW(wc), sO(wc), sG(wc)};
            rates = {qWs, qOs, qGs};
        else
            rhows = [f.rhoOS, f.rhoGS];
            bw    = {bO(wc), bG(wc)};
            mw    = {mobO(wc), mobG(wc)};
            sat   = {sO(wc), sG(wc)};
            rates = {qOs, qGs};
        end
        
        L_ix = woffset + 1;
        V_ix = woffset + 2;
        
        [cqs, weqs, ctrleq, wc, state.wellSol]  = wm.computeWellFlux(model, W, wellSol, ...
            bhp, rates, pw, rhows, bw, mw, sat, {},...
            'nonlinearIteration', opt.iteration);
        
        eqs(offset + (1:(2 + woffset))) = weqs;
        eqs{offset + (3 + woffset)} = ctrleq;
        
        injO = double(cqs{L_ix}) > 0;
        injG = double(cqs{V_ix}) > 0;
        
        if model.water
            % Water
            eqs{1}(wc) = eqs{1}(wc) - fluid.rhoWS.*cqs{1};
            names(offset + (1:4)) = {'waterWells', 'oilWells', 'gasWells', 'closureWells'};
            types(offset + (1:4)) = {'perf', 'perf', 'perf', 'well'};
        else
            names(offset + (1:3)) = {'oilWells', 'gasWells', 'closureWells'};
            types(offset + (1:3)) = {'perf', 'perf', 'well'};
        end
        
        srcTot = 0;
        compSrc = zeros(numel(wc), ncomp);
        for i = 1:ncomp
            ix = i + woffset;
            % Mixture of production and injection. Production uses cell
            % values for components, injectors use whatever was prescribed
            % to the well.
            src =       (fluid.rhoOS.*cqs{L_ix}.*injO + fluid.rhoGS.*cqs{V_ix}.*injG).*w_comp(wm.perf2well, i)...
                       + fluid.rhoOS.*~injO.*x_comp{i}.*cqs{L_ix} + fluid.rhoGS.*~injG.*y_comp{i}.*cqs{V_ix};
            
            eqs{ix}(wc) = eqs{ix}(wc) - src;

            compSrc(:, i) = double(src);
            srcTot = srcTot + double(src);
        end
        fluxt = double(srcTot);
        for i = 1:numel(W)
            wp = wm.perf2well == i;
            state.wellSol(i).flux = fluxt(wp);
            state.wellSol(i).components = compSrc(i, :);
        end
    end
end

if model.water
    wscale = dt./(s.pv*mean(double(rhoW)));
    eqs{1} = eqs{1}.*wscale;
end
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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

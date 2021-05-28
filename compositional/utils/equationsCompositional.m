function [problem, state] = equationsCompositional(state0, state, model, dt, drivingForces, varargin)
% Overall composition fully-implicit equations

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

opt = struct('Verbose',     mrstVerbose,...
            'reverseMode', false,...
            'resOnly',     false,...
            'pressure',    false, ...
            'propsPressure', [], ...
            'staticWells',  false, ...
            'iteration',   -1);

opt = merge_options(opt, varargin{:});

% Shorter names for some commonly used parts of the model and forces.
s = model.operators;
f = model.fluid;
W = drivingForces.W;

fluid = model.fluid;
cnames = model.EOSModel.getComponentNames();

% Properties at current timestep
[p, sW, z, temp, wellSol] = model.getProps(state, ...
    'pressure', 'water', 'z', 'T', 'wellSol');

[p0, sW0, sO0, sG0, z0, temp0, wellSol0] = model.getProps(state0, ...
    'pressure', 'water', 'oil', 'gas', 'z', 'T', 'wellSol');
z = ensureMinimumFraction(z);
z0 = ensureMinimumFraction(z0);
z = expandMatrixToCell(z);
z0 = expandMatrixToCell(z0);

[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

ncomp = numel(z);
cnames = model.EOSModel.getComponentNames();

nwellvar = sum(cellfun(@numel, wellVars));
if model.water
    [p, z{1:ncomp-1}, sW, wellVars{:}] = model.AutoDiffBackend.initVariablesAD(p, z{1:ncomp-1}, sW, wellVars{:});
    primaryVars = {'pressure', cnames{1:end-1}, 'sW', wellVarNames{:}};
else
    [p, z{1:ncomp-1}, wellVars{:}] = model.AutoDiffBackend.initVariablesAD(p, z{1:ncomp-1}, wellVars{:});
    primaryVars = {'pressure', cnames{1:end-1}, wellVarNames{:}};
end
% Property pressure different from flow potential
p_prop = opt.propsPressure;
otherPropPressure = ~isempty(p_prop);
if ~otherPropPressure
    p_prop = p;
end


z{end} = 1;
for i = 1:(ncomp-1)
    z{end} = z{end} - z{i};
end

[xM,  yM,  sO,  sG,  rhoO,  rhoG, muO, muG] = model.computeTwoPhaseFlowProps(state, p_prop, temp, z);
[xM0, yM0, ~, ~, rhoO0, rhoG0] = model.computeTwoPhaseFlowProps(state0, p0, temp0, z0);

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p_prop, p0);

if model.water
    sO = (1-sW).*sO;
    sG = (1-sW).*sG;
    sat = {sW, sO, sG};
    
    [krW, krO, krG] = model.evaluateRelPerm(sat);
    krW = mobMult.*krW;
else
    sat = {sO, sG};
    [krO, krG] = model.evaluateRelPerm(sat);
end

krO = mobMult.*krO;
krG = mobMult.*krG;

% Compute transmissibility
T = s.T.*transMult;

% Gravity gradient per face
gdz = model.getGravityGradient();

% Oil flux
rhoOf  = s.faceAvg(sO.*rhoO)./max(s.faceAvg(sO), 1e-8);
mobO   = krO./muO;
dpO    = s.Grad(p) - rhoOf.*gdz;
upco  = (value(dpO)<=0);
vO = -s.faceUpstr(upco, mobO).*T.*dpO;

% Gas flux
rhoGf  = s.faceAvg(sG.*rhoG)./max(s.faceAvg(sG), 1e-8);
mobG   = krG./muG;
if isfield(fluid, 'pcOG')
    pcOG  = fluid.pcOG(sG);
    pG = p + pcOG;
else
    pG = p;
end
dpG    = s.Grad(pG) - rhoGf.*gdz;
upcg  = (value(dpG)<=0);
vG = -s.faceUpstr(upcg, mobG).*T.*dpG;

rOvO = s.faceUpstr(upco, rhoO).*vO;
rGvG = s.faceUpstr(upcg, rhoG).*vG;


bO = rhoO./fluid.rhoOS;
bG = rhoG./fluid.rhoGS;
% EQUATIONS -----------------------------------------------------------
% water equation + n component equations
[eqs, fluxes, types, names] = deal(cell(1, ncomp + model.water));

woffset = model.water;
if model.water
    % Water flux
    muW = f.muW(p_prop);
    bW     = fluid.bW(p_prop);
    rhoW   = bW.*fluid.rhoWS;
    bW0 = fluid.bW(p0);

    rhoWf  = s.faceAvg(rhoW);
    mobW   = krW./muW;
    if isfield(fluid, 'pcOW')
        pcOW  = fluid.pcOW(sW);
        pW = p - pcOW;
    else
        pW = p;
    end
    dpW    = s.Grad(pW) - rhoWf.*gdz;
    upcw  = (value(dpW)<=0);
    vW = -s.faceUpstr(upcw, mobW).*T.*dpW;
    rWvW = s.faceUpstr(upcw, bW).*vW;
    
    eqs{1} = (s.pv/dt).*( bW.*pvMult.*sW - bW0.*pvMult0.*sW0);
    names{1} = 'water';
    types{1} = 'cell';
    fluxes{1} = rWvW;
else
    [vW, mobW, upcw, bW, rhoW] = deal([]);
end

if model.outputFluxes
    state = model.storeFluxes(state, vW, vO, vG);
end
if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, bG);
    state = model.storeMobilities(state, mobW, mobO, mobG);
    state = model.storeUpstreamIndices(state, upcw, upco, upcg);
    state = model.storeDensities(state, rhoW, rhoO, rhoG);
end

if woffset
end

for i = 1:ncomp
    names{i+woffset} = cnames{i};
    types{i+woffset} = 'cell';

    eqs{i+woffset} = (s.pv/dt).*( ...
                    rhoO.*pvMult.*sO.*xM{i} - rhoO0.*pvMult0.*sO0.*xM0{i} + ...
                    rhoG.*pvMult.*sG.*yM{i} - rhoG0.*pvMult0.*sG0.*yM0{i});
    fluxes{i+woffset} = rOvO.*s.faceUpstr(upco, xM{i}) + rGvG.*s.faceUpstr(upcg, yM{i});
end


% Finally, add in and setup well equations
if model.water
    rho = {rhoW, rhoO, rhoG};
    mob = {mobW, mobO, mobG};
    pressures = {pW, p, pG};
else
    rho = {rhoO, rhoG};
    mob = {mobO, mobG};
    pressures = {p, pG};
end
comps = cellfun(@(x, y) {x, y}, xM, yM, 'UniformOutput', false);


[eqs, state] = model.addBoundaryConditionsAndSources(eqs, names, types, state, ...
                                                 pressures, sat, mob, rho, ...
                                                 {}, comps, ...
                                                 drivingForces);

if opt.staticWells
    compSrc = vertcat(wellSol.components);
    for i = 1:ncomp
        wc = vertcat(W.cells);
        eqs{i+woffset}(wc) = eqs{i+woffset}(wc) - compSrc(:, i);
    end
else
    [eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, ...
                                                     types, wellSol0, wellSol, ...
                                                     wellVars, ...
                                                     wellMap, p, mob, rho, ...
                                                     {}, comps, ...
                                                     dt, opt);
end

conservationindices = 1:(ncomp+woffset);
acc = eqs(conservationindices);
for i = conservationindices
    eqs{i} = s.AccDiv(eqs{i}, fluxes{i});
    if model.water && i > woffset
        pureWater = value(sW) == 1;
        if any(pureWater)
            % Cells with pure water should just retain their composition to
            % avoid singular systems
            eqs{i}(pureWater) = eqs{i}(pureWater) + ...
                            1e-3*(z{i}(pureWater) - value(z{i}(pureWater)));
        end
    end
end

if opt.pressure
    wat = model.water;
    if opt.resOnly
        weights = cell(1, ncomp + wat);
        [weights{:}] = deal(1);
    else
        state = model.storeDensities(state, rhoW, rhoO, rhoG);
        ncell = model.G.cells.num;
        ndof = ncell*(ncomp+wat);
        wellVarIndices = (ndof+1):(ndof+nwellvar);
        
        [w, dwdp] = getPartialVolumes(model, state, acc, ...
            'iteration',                  opt.iteration, ...
            'wellVarIndices',             wellVarIndices, ...
            'singlePhaseStrategy',        model.singlePhaseStrategy, ...
            'twoPhaseStrategy',           model.twoPhaseStrategy, ...
            'singlePhaseDifferentiation', model.singlePhaseDifferentiation, ...
            'twoPhaseDifferentiation',    model.twoPhaseDifferentiation);

        weights = cell(ncomp+wat, 1);
        for i = 1:(ncomp+wat)
            wi = w(:, i);
            if any(dwdp)
                Wp = double2ADI(wi, p);
                Wp.jac{1} = sparse(1:ncell, 1:ncell, dwdp(:, i), ncell, ncell);
            else
                Wp = wi;
            end
            weights{i} = Wp;
        end
        state.w = w;
        state.w_p = state.pressure;
    end

    
    peq = 0;
    for i = 1:(ncomp+wat)
        peq = peq + weights{i}.*eqs{i};
    end
    active = false(numel(eqs), 1);
    active(1) = true;
    active(ncomp+wat+1:end) = true;

    eqs{1} = peq;
    
    eqs = eqs(active);
    for i = 1:numel(eqs)
        eqs{i}.jac = eqs{i}.jac(active);
    end
    names{1} = 'pressure';
    primaryVars = primaryVars(active);
    names = names(active);
    types = types(active);
end

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
problem.iterationNo = opt.iteration;
end

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

function [problem, state] = equationsNaturalVariables(state0, state, model, dt, drivingForces, varargin)
% Equations for natural variables formulation

%{
Copyright 2009-2017 SINTEF Digital, Mathematics & Cybernetics.

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
            'staticWells',  false, ...
            'propsPressure', [], ...
            'reduceToPressure', false, ...
            'iteration',   -1);

opt = merge_options(opt, varargin{:});
% Shorter names for some commonly used parts of the model and forces.
s = model.operators;
f = model.fluid;
W = drivingForces.W;

fluid = model.fluid;
compFluid = model.EOSModel.fluid;
% Properties at current timestep
[p, sW, sO, sG, x, y, z, temp, wellSol] = model.getProps(state, ...
    'pressure', 'water', 'so', 'sg', 'x', 'y', 'z', 'T', 'wellSol');
z = expandMatrixToCell(z);
[p0, sW0, sO0, sG0, x0, y0, temp0, wellSol0] = model.getProps(state0, ...
    'pressure', 'water', 'so', 'sg', 'x', 'y', 'T', 'wellSol');
[pureLiquid, pureVapor, twoPhase] = model.getFlag(state);

multiPhase = twoPhase;
freeOil = twoPhase;
freeGas = twoPhase;

x = expandMatrixToCell(x);
y = expandMatrixToCell(y);

[wellvars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);
if opt.staticWells
    wellvars = {[]};
end
ncomp = compFluid.getNumberOfComponents();
[xnames, ynames, cnames] = deal(model.EOSModel.fluid.names);
for i = 1:ncomp
    xnames{i} = ['v_', cnames{i}];
    ynames{i} = ['w_', cnames{i}];
end

twoPhaseIx = find(twoPhase);

wtmp = ones(nnz(twoPhase), 1);

w = cell(ncomp, 1);
[w{:}] = deal(wtmp);

nc = model.G.cells.num;
for i = 1:(ncomp-1)
    w{i} = y{i}(twoPhase);
end

so = sO(freeOil);

nwellvar = sum(cellfun(@numel, wellvars));
nwelleqs = numel(wellvars);


if opt.resOnly
    initfn = @deal;
else
    initfn = @(varargin) model.AutoDiffBackend.initVariablesAD(varargin{:});
end
sg = sG(freeGas);

if model.water
    [p, x{1:ncomp-1}, sW, wellvars{:}, so, w{1:ncomp-1}, sg] = initfn(...
     p, x{1:ncomp-1}, sW, wellvars{:}, so, w{1:ncomp-1}, sg);
    primaryVars = {'pressure', xnames{1:end-1}, 'satw', wellVarNames{:}, 'sato', ynames{1:end-1}, 'satg'};
else
    [p, x{1:ncomp-1}, wellvars{:}, so, w{1:ncomp-1}, sg] = initfn(...
     p, x{1:ncomp-1}, wellvars{:}, so, w{1:ncomp-1}, sg);
    primaryVars = {'pressure', xnames{1:end-1}, wellVarNames{:}, 'sato', ynames{1:end-1}, 'satg'};
    sW = zeros(model.G.cells.num, 1);
end
sample = p;

sO = model.AutoDiffBackend.convertToAD(sO, sample);
sO(freeOil) = so;
sG = model.AutoDiffBackend.convertToAD(sG, sample);
sG(freeGas) = sg;


if model.water
    if any(pureVapor)
        sG(pureVapor) = 1 - sW(pureVapor);
    end
    
    if any(pureLiquid)
        sO(pureLiquid) = 1 - sW(pureLiquid);
    end
end

if isempty(twoPhaseIx) || opt.resOnly
    reorder = [];
else
    n2ph = nnz(twoPhase);
    nVars = sum(p.getNumVars());
    reorder = 1:nVars;
    start = nc + twoPhaseIx;
    stop = nc*(ncomp+model.water) + nwellvar + (1:n2ph);
    reorder(start) = stop;
    reorder(stop) = start;
end


% Property pressure different from flow potential
p_prop = opt.propsPressure;
otherPropPressure = ~isempty(p_prop);
if ~otherPropPressure
    p_prop = p;
end


cellJacMap = cell(numel(primaryVars), 1);

offset = ncomp + model.water + numel(wellvars);
if any(twoPhase) && ~all(twoPhase)
    for i = 1:(ncomp+1)
        cellJacMap{i + offset} = twoPhaseIx;
    end
end
x{end} = ones(model.G.cells.num, 1);
w{end} = ones(nnz(twoPhase), 1);

for i = 1:ncomp-1
    x{end} = x{end}-x{i};
    if any(twoPhase)
        w{end} = w{end}-w{i};
    end
end

for i = 1:ncomp
    y{i} = ~pureLiquid.*x{i} + double(x{i}).*pureLiquid;
    if any(twoPhase)
        if ~opt.resOnly
            assert(isa(y{i}, 'ADI'));
        end
        y{i}(twoPhase) = w{i};
    end
    x{i}(pureVapor) = double(x{i}(pureVapor));
end


% Compute properties and fugacity
[xM,  yM,  rhoO,  rhoG,  muO,  muG, f_L, f_V, xM0, yM0, rhoO0, rhoG0] = ...
                  model.getTimestepPropertiesEoS(state, state0, p_prop, temp, x, y, z, sO, sG, cellJacMap);

if model.water
    sat = {sW, sO, sG};
else
    sat = {sO, sG};
end


if model.water
    [krW, krO, krG] = model.evaluateRelPerm(sat);
else
    [krO, krG] = model.evaluateRelPerm(sat);
end

% Compute transmissibility
T = s.T;

% Gravity gradient per face
gdz = model.getGravityGradient();

% Oil flux
rhoOf  = s.faceAvg(sat{1+model.water}.*rhoO)./max(s.faceAvg(sat{1+model.water}), 1e-8);
mobO   = krO./muO;
dpO    = s.Grad(p) - rhoOf.*gdz;
upco  = (double(dpO)<=0);
vO = -s.faceUpstr(upco, mobO).*T.*dpO;

% Gas flux
rhoGf  = s.faceAvg(sat{2+model.water}.*rhoG)./max(s.faceAvg(sat{2+model.water}), 1e-8);
mobG   = krG./muG;
dpG    = s.Grad(p) - rhoGf.*gdz;
upcg  = (double(dpG)<=0);
vG = -s.faceUpstr(upcg, mobG).*T.*dpG;

rOvO = s.faceUpstr(upco, rhoO).*vO;
rGvG = s.faceUpstr(upcg, rhoG).*vG;

pv = model.operators.pv;
pv0 = pv;
if isfield(fluid, 'pvMultR')
    pv = pv.*fluid.pvMultR(p_prop);
    pv0 = pv0.*fluid.pvMultR(p0);
end


% EQUATIONS -----------------------------------------------------------
if model.water
    % Water flux
    muW = f.muW(p_prop);
    bW     = fluid.bW(p_prop);
    rhoW   = bW.*fluid.rhoWS;
    rhoW0 = fluid.bW(p0).*fluid.rhoWS;

    rhoWf  = s.faceAvg(rhoW);
    mobW   = krW./muW;
    dpW    = s.Grad(p) - rhoWf.*gdz;
    upcw  = (double(dpW)<=0);
    vW = -s.faceUpstr(upcw, mobW).*T.*dpW;
    rWvW = s.faceUpstr(upcw, rhoW).*vW;
    water = (1/dt).*(pv.*rhoW.*sW - pv0.*rhoW0.*sW0);
    divWater = s.Div(rWvW);
else
    [vW, mobW, upcw, bW, rhoW] = deal([]);
end

if model.outputFluxes
    state = model.storeFluxes(state, vW, vO, vG);
end

if model.extraStateOutput
    bO = rhoO./fluid.rhoOS;
    bG = rhoG./fluid.rhoGS;
    state = model.storebfactors(state, bW, bO, bG);
    state = model.storeMobilities(state, mobW, mobO, mobG);
    state = model.storeUpstreamIndices(state, upcw, upco, upcg);
    state = model.storeDensities(state, rhoW, rhoO, rhoG);
end

state = model.storeDensities(state, rhoW, rhoO, rhoG);
% water equation + n component equations
[eqs, types, names] = deal(cell(1, ncomp + model.water));


if opt.reduceToPressure
    C = cell(2*ncomp + model.water, 1);
end

divTerm = cell(ncomp, 1);
for i = 1:ncomp
    names{i} = compFluid.names{i};
    types{i} = 'cell';
    eqs{i} = (1/dt).*( ...
                    pv.*rhoO.*sO.*xM{i} - pv0.*rhoO0.*sO0.*xM0{i} + ...
                    pv.*rhoG.*sG.*yM{i} - pv0.*rhoG0.*sG0.*yM0{i});
    
    divTerm{i} = s.Div(rOvO.*s.faceUpstr(upco, xM{i}) + rGvG.*s.faceUpstr(upcg, yM{i}));
    if opt.reduceToPressure
        C{i} = eqs{i};
    end
    if model.water
        pureWater = double(sW) == 1;
        if any(pureWater)
            % Cells with pure water should just retain their composition to
            % avoid singular systems
            eqs{i}(pureWater) = eqs{i}(pureWater) + ...
                            1e-3*(x{i}(pureWater) - double(x{i}(pureWater)));
        end
    end
    
end

if model.water
    if ~opt.reduceToPressure
        eqs{ncomp+1} = water.*(dt./(s.pv.*double(rhoW)));
    end
    eqs{ncomp+1} = water;
    names{ncomp+1} = 'water';
    types{ncomp+1} = 'cell';
    C{ncomp+1} = water;
    
    rho = {rhoW, rhoO, rhoG};
    mob = {mobW, mobO, mobG};
    pressures = {p, p, p};
else
    rho = {rhoO, rhoG};
    mob = {mobO, mobG};
    pressures = {p, p};
end
comps = cellfun(@(x, y) {x, y}, xM, yM, 'UniformOutput', false);

[eqs, state] = model.addBoundaryConditionsAndSources(eqs, names, types, state, ...
                                                 pressures, sat, mob, rho, ...
                                                 {}, comps, ...
                                                 drivingForces);

% Finally, add in and setup well equations
if opt.staticWells
    compSrc = vertcat(wellSol.components);
    for i = 1:ncomp
        wc = vertcat(W.cells);
        eqs{i+model.water}(wc) = eqs{i+model.water}(wc) - compSrc(:, i);
    end
    nwelleqs = 0;
else
    wellSol = state.wellSol;
    [eqs, names, types, state.wellSol, src] = model.insertWellEquations(eqs, names, ...
                                                     types, wellSol0, wellSol, ...
                                                     wellvars, ...
                                                     wellMap, p, mob, rho, ...
                                                     {}, comps, ...
                                                     dt, opt); %#ok
end

eq_offset = nwelleqs + ncomp + model.water;
for i = 1:ncomp
    eqs{i} = eqs{i} + divTerm{i};
    ix = i + eq_offset;
    names{ix}= ['f_', compFluid.names{i}];
    types{ix} = 'fugacity';
    eqs{ix} = (f_L{i}(twoPhase) - f_V{i}(twoPhase))/barsa;
    
    if opt.reduceToPressure
        C{ix - nwelleqs} = eqs{ix};
    end
end

if model.water
    eqs{ncomp+1} = eqs{ncomp+1} + divWater;
end

cloix = eq_offset + ncomp + 1;

if any(multiPhase)
    eqs{cloix} = sW(multiPhase) + sO(multiPhase) + sG(multiPhase) - 1;
else
    eqs{cloix} = [];
end

types{cloix} = 'saturation';
names{cloix} = 'volclosure';

if opt.reduceToPressure
    C{cloix - nwelleqs} = eqs{cloix};
    problem = PressureReducedLinearSystem(eqs, types, names, primaryVars, state, dt);
    problem.accumulationTerms = C;
    problem.model = model;
    problem.wellVarIndices = nc*(ncomp+model.water) + (1:nwellvar);
    problem.nwellvar = nwelleqs;
    problem.wellvars = wellvars;
    problem.wellvarNames = wellVarNames;
else
    scale = (dt./s.pv)./mean(double(sO0).*double(rhoO0) + double(sG0).*double(rhoG0));
    for i = 1:ncomp
        eqs{i} = eqs{i}.*scale;
    end
    if model.reduceLinearSystem
        problem = ReducedLinearizedSystem(eqs, types, names, primaryVars, state, dt);
    else
        problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
    end
end

if isa(problem, 'ReducedLinearizedSystem')
    % problem.keepNum = model.G.cells.num*(ncomp+model.water);
    problem.keepNum = nc*(ncomp+model.water) + nwellvar;
    problem.reorder = reorder;
end

problem.iterationNo = opt.iteration;
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

function [problem, state] = equationsNaturalVariablesDP(state0, state, model, dt, drivingForces, varargin)
% Equations for natural variables formulation

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
% Wei
W = drivingForces.W;

fluid = model.fluid;
mixture = model.EOSModel.CompositionalMixture;
% Properties at current timestep
[p, sW, sO, sG, x, y, z, temp, wellSol] = model.getProps(state, ...
    'pressure', 'sw', 'so', 'sg', 'x', 'y', 'z', 'T', 'wellSol');
% Wei
[p_matrix, sW_matrix, sO_matrix, sG_matrix, x_matrix, y_matrix, z_matrix] = model.getProps(state.matrix, ...
    'pressure', 'sw', 'so', 'sg', 'x', 'y', 'z');
z = expandMatrixToCell(z);
z_matrix = ensureMinimumFraction(z_matrix);
[p0, sW0, sO0, sG0, x0, y0, temp0, wellSol0] = model.getProps(state0, ...
    'pressure', 'sw', 'so', 'sg', 'x', 'y', 'T', 'wellSol');
% Wei
[p0_matrix, sW0_matrix, sO0_matrix, sG0_matrix, x0_matrix, y0_matrix] = model.getProps(state0.matrix, ...
    'pressure', 'sw', 'so', 'sg', 'x', 'y');
[pureLiquid, pureVapor, twoPhase] = model.getFlag(state);
% Wei
[pureLiquid_matrix, pureVapor_matrix, twoPhase_matrix] = model.getFlag(state.matrix);
if 1
    stol = 1e-6;
    pureWater = sO + sG < stol;
    sO(~pureVapor & pureWater) = stol;
    sG(~pureLiquid & pureWater) = stol;
    
    [pureLiquid0, pureVapor0, twoPhase0] = model.getFlag(state0);
    pureWater0 = sO0 + sG0 < stol;
    sO0(~pureVapor0 & pureWater0) = stol;
    sG0(~pureLiquid0 & pureWater0) = stol;
    % Wei
    pureWater_matrix = sO_matrix + sG_matrix < stol;
    sO_matrix(~pureVapor_matrix & pureWater_matrix) = stol;
    sG_matrix(~pureLiquid_matrix & pureWater_matrix) = stol;
    
    [pureLiquid0_matrix, pureVapor0_matrix, twoPhase0_matrix] = model.getFlag(state0.matrix);
    pureWater0_matrix = sO0_matrix + sG0_matrix < stol;
    sO0_matrix(~pureVapor0_matrix & pureWater0_matrix) = stol;
    sG0_matrix(~pureLiquid0_matrix & pureWater0_matrix) = stol;
end

multiPhase = twoPhase;
freeOil = twoPhase;
freeGas = twoPhase;
% Wei
multiPhase_matrix = twoPhase_matrix;
freeOil_matrix = twoPhase_matrix;
freeGas_matrix = twoPhase_matrix;

z_tol = model.EOSModel.minimumComposition;
state0.x = ensureMinimumFraction(state0.x, z_tol);
state0.y = ensureMinimumFraction(state0.y, z_tol);
state0.components = ensureMinimumFraction(state0.components, z_tol);
% Wei
state0.matrix.x = ensureMinimumFraction(state0.matrix.x, z_tol);
state0.matrix.y = ensureMinimumFraction(state0.matrix.y, z_tol);
state0.matrix.components = ensureMinimumFraction(state0.matrix.components, z_tol);

x = ensureMinimumFraction(x, z_tol);
y = ensureMinimumFraction(y, z_tol);
x = expandMatrixToCell(x);
y = expandMatrixToCell(y);
% Wei
x_matrix = ensureMinimumFraction(x_matrix, z_tol);
y_matrix = ensureMinimumFraction(y_matrix, z_tol);
x_matrix = expandMatrixToCell(x_matrix);
y_matrix = expandMatrixToCell(y_matrix);

[wellvars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);
if opt.staticWells
    wellvars = {[]};
end
ncomp = mixture.getNumberOfComponents();
% Wei
[xnames, ynames, cnames, xnames_matrix, ynames_matrix, cnames_matrix] = deal(mixture.names);
for i = 1:ncomp
    xnames{i} = ['v_', cnames{i}];
    ynames{i} = ['w_', cnames{i}];
    % Wei
    xnames_matrix{i} = ['v_', cnames_matrix{i},'_matrix'];
    ynames_matrix{i} = ['w_', cnames_matrix{i},'_matrix'];
end

twoPhaseIx = find(twoPhase);
% Wei
twoPhaseIx_matrix = find(twoPhase_matrix);

wtmp = ones(nnz(twoPhase), 1);
% Wei
wtmp_matrix = ones(nnz(twoPhase_matrix), 1);
% Wei
w = cell(ncomp, 1);
w_matrix = cell(ncomp, 1);
[w{:}] = deal(wtmp);
[w_matrix{:}] = deal(wtmp_matrix);

nc = model.G.cells.num;
for i = 1:(ncomp-1)
    w{i} = y{i}(twoPhase);
    w_matrix{i} = y_matrix{i}(twoPhase_matrix);
end

so = sO(freeOil);
% Wei
so_matrix = sO_matrix(freeOil_matrix);

nwellvar = sum(cellfun(@numel, wellvars));
nwelleqs = numel(wellvars);

if opt.resOnly
    initfn = @deal;
else
    initfn = @(varargin) model.AutoDiffBackend.initVariablesAD(varargin{:});
end
sg = sG(freeGas);
sg_matrix = sG_matrix(freeGas_matrix);
% Wei
if model.water
    [p, x{1:ncomp-1}, p_matrix, x_matrix{1:ncomp-1}, sW, sW_matrix, wellvars{:}, so, w{1:ncomp-1}, so_matrix, w_matrix{1:ncomp-1}, sg, sg_matrix] = initfn(...
     p, x{1:ncomp-1}, p_matrix, x_matrix{1:ncomp-1}, sW, sW_matrix, wellvars{:}, so, w{1:ncomp-1}, so_matrix, w_matrix{1:ncomp-1}, sg, sg_matrix);
    primaryVars = {'pressure', xnames{1:end-1}, 'pressure_matrix', xnames_matrix{1:end-1}, 'sw', 'sw_matrix',wellVarNames{:}, 'sl', ynames{1:end-1}, 'sl_matrix', ynames_matrix{1:end-1}, 'sv', 'sv_matrix'};
else
    [p, x{1:ncomp-1}, p_matrix, x_matrix{1:ncomp-1}, wellvars{:}, so, w{1:ncomp-1}, so_matrix, w_matrix{1:ncomp-1}, sg, sg_matrix] = initfn(...
     p, x{1:ncomp-1}, p_matrix, x_matrix{1:ncomp-1}, wellvars{:}, so, w{1:ncomp-1}, so_matrix, w_matrix{1:ncomp-1}, sg, sg_matrix);
    primaryVars = {'pressure', xnames{1:end-1}, 'pressure_matrix', xnames_matrix{1:end-1}, wellVarNames{:}, 'sl', ynames{1:end-1}, 'sl_matrix', ynames_matrix{1:end-1}, 'sv', 'sv_matrix'};
    sW = zeros(model.G.cells.num, 1);
    sW_matrix = zeros(model.G.cells.num, 1);
end
sample = p;
sample_matrix = p_matrix;

sO = model.AutoDiffBackend.convertToAD(sO, sample);
sO(freeOil) = so;
sG = model.AutoDiffBackend.convertToAD(sG, sample);
sG(freeGas) = sg;
% Wei
sO_matrix = model.AutoDiffBackend.convertToAD(sO_matrix, sample_matrix);
sO_matrix(freeOil_matrix) = so_matrix;
sG_matrix = model.AutoDiffBackend.convertToAD(sG_matrix, sample_matrix);
sG_matrix(freeGas_matrix) = sg_matrix;

[sO, sG] = setMinimums(model, state, sW, sO, sG, pureVapor, pureLiquid);
[sO0, sG0] = setMinimums(model, state0, sW0, sO0, sG0, pureVapor0, pureLiquid0);
% Wei
[sO_matrix, sG_matrix] = setMinimums(model, state.matrix, sW_matrix, sO_matrix, sG_matrix, pureVapor_matrix, pureLiquid_matrix);
[sO0_matrix, sG0_matrix] = setMinimums(model, state0.matrix, sW0_matrix, sO0_matrix, sG0_matrix, pureVapor0_matrix, pureLiquid0_matrix);

if isempty(twoPhaseIx) || opt.resOnly
    reorder = [];
else
    n2ph = nnz(twoPhase);
    nVars = sum(p.getNumVars());
    reorder = 1:nVars;
    start = nc + twoPhaseIx;
    stop = nc*(2*ncomp+2*model.water) + nwellvar + (1:n2ph);
    reorder(start) = stop;
    reorder(stop) = start;
end
% Wei
if isempty(twoPhaseIx_matrix) || opt.resOnly
    reorder = [];
else
    n2ph_matrix = nnz(twoPhase_matrix);
    if isempty(reorder)
        nVars = sum(p.getNumVars());
        reorder = 1:nVars;
    end
    start = nc*(ncomp+1) + twoPhaseIx_matrix;
    stop = nc*(2*ncomp+2*model.water) + nwellvar + ncomp*n2ph_matrix + (1:n2ph_matrix);
    reorder(start) = stop;
    reorder(stop) = start;
end

% Property pressure different from flow potential
p_prop = opt.propsPressure;
otherPropPressure = ~isempty(p_prop);
if ~otherPropPressure
    p_prop = p;
end
% Wei
p_prop_matrix = p_matrix;

cellJacMap = cell(numel(primaryVars), 1);
cellJacMap_matrix = cellJacMap;

offset = ncomp + model.water + numel(wellvars);
if any(twoPhase) && ~all(twoPhase)
    for i = 1:(ncomp+1)
        cellJacMap{i + offset} = twoPhaseIx;
    end
end
if any(twoPhase_matrix) && ~all(twoPhase_matrix)
    for i = 1:(ncomp+1)
        cellJacMap_matrix{i + offset} = twoPhaseIx_matrix;
    end
end
x{end} = ones(model.G.cells.num, 1);
w{end} = ones(nnz(twoPhase), 1);
% Wei
x_matrix{end} = ones(model.G.cells.num, 1);
w_matrix{end} = ones(nnz(twoPhase_matrix), 1);
for i = 1:ncomp-1
    x{end} = x{end}-x{i};
    x_matrix{end} = x_matrix{end}-x_matrix{i};
    if any(twoPhase)
        w{end} = w{end}-w{i};
        w_matrix{end} = w_matrix{end}-w_matrix{i};
    end
end

for i = 1:ncomp
    y{i} = ~pureLiquid.*x{i} + value(x{i}).*pureLiquid;
    y_matrix{i} = ~pureLiquid_matrix.*x_matrix{i} + value(x_matrix{i}).*pureLiquid_matrix;
    if any(twoPhase)
        if ~opt.resOnly
            assert(isa(y{i}, 'ADI'));
        end
        y{i}(twoPhase) = w{i};
    end
    x{i}(pureVapor) = value(x{i}(pureVapor));
    % Wei
    if any(twoPhase_matrix)
        if ~opt.resOnly
            assert(isa(y_matrix{i}, 'ADI'));
        end
        y_matrix{i}(twoPhase_matrix) = w_matrix{i};
    end
    x_matrix{i}(pureVapor_matrix) = value(x_matrix{i}(pureVapor_matrix));
end


% Compute properties and fugacity
[xM,  yM,  rhoO,  rhoG,  muO,  muG, f_L, f_V, xM0, yM0, rhoO0, rhoG0] = ...
                  model.getTimestepPropertiesEoS(state, state0, p_prop, temp, x, y, z, sO, sG, cellJacMap);
% Wei
[xM_matrix,  yM_matrix,  rhoO_matrix,  rhoG_matrix,  muO_matrix,  muG_matrix, f_L_matrix, f_V_matrix, xM0_matrix, yM0_matrix, rhoO0_matrix, rhoG0_matrix] = ...
                  model.getTimestepPropertiesEoS(state.matrix, state0.matrix, p_prop_matrix, temp, x_matrix, y_matrix, z_matrix, sO_matrix, sG_matrix, cellJacMap_matrix);
if model.water
    sat = {sW, sO, sG};
    sat_matrix = {sW_matrix, sO_matrix, sG_matrix};

    [krW, krO, krG] = model.evaluateRelPerm(sat);
else
    sat = {sO, sG};
    sat_matrix = {sO_matrix, sG_matrix};
    if ~isfield(model.fluid, 'relPermScal')
        [krO, krG] = model.evaluateRelPerm(sat,'medium','fracture');
    else
        % Wei edit Relative Permeability Scaling, but some function is missing
        kr = eCPARelativePermeabilityScaling(model, state);
        [krO, krG] = deal(kr{:});
    end
end

% Compute transmissibility
T = s.T;

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
    % Wei
    pcOG_matrix  = fluid.pcOG(sG_matrix);
    pG_matrix = p_matrix + pcOG_matrix;
else
    pG = p;
    % Wei
    pG_matrix = p_matrix;
end
dpG    = s.Grad(pG) - rhoGf.*gdz;

upcg  = (value(dpG)<=0);
vG = -s.faceUpstr(upcg, mobG).*T.*dpG;

rOvO = s.faceUpstr(upco, rhoO).*vO;
rGvG = s.faceUpstr(upcg, rhoG).*vG;

pv = model.operators.pv;
pv0 = pv;
% Wei
pv_matrix = model.operators.pv_matrix;
pv0_matrix = pv_matrix;
if isfield(fluid, 'pvMultR')
    pv = pv.*fluid.pvMultR(p_prop);
    pv0 = pv0.*fluid.pvMultR(p0);
    % Wei
    pv_matrix = pv_matrix.*fluid.pvMultR(p_prop_matrix);
    pv0_matrix = pv0_matrix.*fluid.pvMultR(p0_matrix);
end

% EQUATIONS -----------------------------------------------------------
if model.water
    % Water flux
    if isfield(fluid, 'pcOW')
        pcOW  = fluid.pcOW(sW);
    else
        pcOW = 0;
    end
    pW = p_prop;
    pW0 = p0;
    muW = fluid.muW(pW);
    bW     = fluid.bW(pW);
    rhoW   = bW.*fluid.rhoWS;
    bW0 = fluid.bW(pW0);

    rhoWf  = s.faceAvg(rhoW);
    mobW   = krW./muW;

    dpW    = s.Grad(p - pcOW) - rhoWf.*gdz;
    upcw  = (value(dpW)<=0);
    vW = -s.faceUpstr(upcw, mobW).*T.*dpW;
    rWvW = s.faceUpstr(upcw, bW).*vW;
    water = (1/dt).*(pv.*bW.*sW - pv0.*bW0.*sW0);
else
    [vW, mobW, upcw, bW, rhoW, rhoW_matrix] = deal([]);
end
%% Wei: Transfer
transfer_model = model.transfer_model_object;
sigma = transfer_model.shape_factor_object.calculate_shape_factor(model);
%% Wei: Fluid heights
%vertical matrix block height
lz=transfer_model.shape_factor_object.block_dimension(:,3);

hG = sG .* lz;
hG_matrix = sG_matrix .* lz;
hO = sO .* lz;
hO_matrix = sO_matrix .* lz;

% Wei: here the matrix densities have been used
delta_rhoGO_matrix = (rhoG_matrix-rhoO_matrix);

%Gravity
g = norm(gravity);

% Introduce pseudo-pressures to account for the pressure
% gradient due to gravity, see Eclipse Technical Description
psi_O = p + delta_rhoGO_matrix.*g.*hO;
psi_G = pG + delta_rhoGO_matrix.*g.*hG;
psi_O_matrix = p_matrix + delta_rhoGO_matrix.*g.*hO_matrix;
psi_G_matrix = pG_matrix + delta_rhoGO_matrix.*g.*hG_matrix;

dpsiG = (value(psi_G_matrix-psi_G)<=0);
dpsiO = (value(psi_O_matrix-psi_O)<=0);

[krGt, krOt, TG, TO, Tdg] = deal(cell(1, ncomp));
for i = 1:ncomp
% Upstream weighing of mobilities
krGt{i} = krG.*yM{i}.*rhoG./muG.*dpsiG + krG.*yM_matrix{i}.*rhoG_matrix./muG_matrix.*(~dpsiG);
krOt{i} = krO.*xM{i}.*rhoO./muO.*dpsiO + krO.*xM_matrix{i}.*rhoO_matrix./muO_matrix.*(~dpsiO);

%% Compute Transfer
%(units 1/T)
vb = model.G.cells.volumes;
TG{i} = vb.*(sigma.*krGt{i}).*(psi_G - psi_G_matrix);
TO{i} = vb.*(sigma.*krOt{i}).*(psi_O - psi_O_matrix);
end

% Wei: Matrix-fracture diffusion
if isfield(model.rock, 'Mechanisms') && isfield(model.rock.Mechanisms,'diffusion')
    kx = model.rock_matrix.perm(:,1);
    try
        ky = model.rock_matrix.perm(:,2);
        kz = model.rock_matrix.perm(:,3);
    catch
        ky = model.rock_matrix.perm(:,1);
        kz = model.rock_matrix.perm(:,1);
    end
    km = (kx.*ky.*kz).^(1/3);

    % dsgG = (value(sG_matrix-sG)<=0);
    % Sg_poro_rho = sG_matrix.*rhoG_matrix.*model.rock_matrix.poro.*dsgG + sG.*rhoG.*model.rock.poro.*(~dsgG);
    % dsgO = (value(sO_matrix-sO)<=0);
    % So_poro_rho = sO_matrix.*rhoO_matrix.*model.rock_matrix.poro.*dsgO + sO.*rhoO.*model.rock.poro.*(~dsgO);

    Sg_poro_rho = sG_matrix.*rhoG_matrix.*model.rock_matrix.poro;
    % So_poro_rho = sO_matrix.*rhoO_matrix.*model.rock_matrix.poro;
    Di = model.rock.Di;
    for i = 1:ncomp        
        Tdg{i} = vb.*sigma./km.*Sg_poro_rho.*Di(i).*(yM{i}-yM_matrix{i});
        % Tdo{i} = vb.*sigma./km.*So_poro_rho.*Di(i).*(xM{i}-xM_matrix{i});
    end
else
    for i = 1:ncomp
        Tdg{i} = 0;
        % Tdo{i} = 0;
    end
end

% Wei
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
state.matrix = model.storeDensities(state.matrix, rhoW_matrix, rhoO_matrix, rhoG_matrix);
% water equation + n component equations
[eqs, types, names] = deal(cell(1, 2*ncomp + model.water));

if opt.reduceToPressure
    C = cell(2*ncomp + model.water, 1);
end

fluxes = cell(2*ncomp, 1);
% Wei
compFlux = zeros(model.G.faces.num, ncomp);

%Wei edit: modified the accumulation term in the mass balance equation
    %to include the contributions of the adsorption in shale formation.
if isfield(model.rock, 'Mechanisms') && isfield(model.rock.Mechanisms,'sorption')
    sumIso = yM_matrix{1}./model.rock.isothermP(1);
    sumIso0 = yM0_matrix{1}./model.rock.isothermP(1);
    for ii=2:numel(yM_matrix)
        sumIso = sumIso + yM_matrix{ii}./model.rock.isothermP(ii);
        sumIso0 = sumIso0 + yM0_matrix{ii}./model.rock.isothermP(ii);
    end
end
for i = 1:ncomp
    names{i} = mixture.names{i};
    types{i} = 'cell';
    eqs{i} = (1/dt).*( ...
                    pv.*rhoO.*sO.*xM{i} - pv0.*rhoO0.*sO0.*xM0{i} + ...
                    pv.*rhoG.*sG.*yM{i} - pv0.*rhoG0.*sG0.*yM0{i});
    vi = rOvO.*s.faceUpstr(upco, xM{i}) + rGvG.*s.faceUpstr(upcg, yM{i});
    fluxes{i} = vi;
    compFlux(model.operators.internalConn,i) = value(vi);

    % Wei
    names{ncomp+i} = [mixture.names{i},'_matrix'];
    types{ncomp+i} = 'cell';
    if isfield(model.rock, 'Mechanisms') && isfield(model.rock.Mechanisms,'sorption')
        eqs{ncomp+i} = (1/dt).*( ...
            pv_matrix.*rhoO_matrix.*sO_matrix.*xM_matrix{i} - pv0_matrix.*rhoO0_matrix.*sO0_matrix.*xM0_matrix{i} + ...
            pv_matrix.*rhoG_matrix.*sG_matrix.*yM_matrix{i} - pv0_matrix.*rhoG0_matrix.*sG0_matrix.*yM0_matrix{i} +...
            (s.gv.*model.rock.isothermRho(i)./model.rock.isothermP(i)).*...
            ((yM_matrix{i}.*p_matrix)./(1+p_matrix.*sumIso) - (yM0_matrix{i}.*p0_matrix)./(1+p0_matrix.*sumIso0)) );
    else
        eqs{ncomp+i} = (1/dt).*( ...
            pv_matrix.*rhoO_matrix.*sO_matrix.*xM_matrix{i} - pv0_matrix.*rhoO0_matrix.*sO0_matrix.*xM0_matrix{i} + ...
            pv_matrix.*rhoG_matrix.*sG_matrix.*yM_matrix{i} - pv0_matrix.*rhoG0_matrix.*sG0_matrix.*yM0_matrix{i});
    end
    fluxes{ncomp+i} = 0.*vi;

    if opt.reduceToPressure
        C{i} = eqs{i};
    end
end
state.componentFluxes = compFlux;

if model.water
    mf = [model.fluid.rhoWS.*value(rWvW), value(rOvO), value(rGvG)];
else
    mf = [value(rOvO), value(rGvG)];
end
state.massFlux = zeros(model.G.faces.num, 2 + model.water);
state.massFlux(model.operators.internalConn, :) = mf;


if model.water
    eqs{ncomp+1} = water;
    names{ncomp+1} = 'water';
    types{ncomp+1} = 'cell';
    C{ncomp+1} = water;
    
    rho = {rhoW, rhoO, rhoG};
    mob = {mobW, mobO, mobG};
    pressures = {pW, p, pG};
else
    rho = {rhoO, rhoG};
    mob = {mobO, mobG};
    pressures = {p, pG};
    % Wei
    rho_matrix = {rhoO_matrix, rhoG_matrix};
    pressures_matrix = {p_matrix, pG_matrix};
end
comps = cellfun(@(x, y) {x, y}, xM, yM, 'UniformOutput', false);
comps_matrix = cellfun(@(x, y) {x, y}, xM_matrix, yM_matrix, 'UniformOutput', false);

woffset = model.water;
[eqsBC, state, bsrc] = model.addBoundaryConditionsAndSources(eqs(1:ncomp+woffset), names(1:ncomp+woffset), types(1:ncomp+woffset), state, ...
                                                 pressures, sat, mob, rho, ...
                                                 {}, comps, ...
                                                 drivingForces);
eqs(1:ncomp+woffset) = eqsBC;

if ~isempty(drivingForces.bc)
    % Compute the contribution from BC while avoiding the capillary jump

    [eqsmat, state.matrix] = addBoundaryConditionsAndSources(model, ...
        eqs(1+ncomp+woffset:2*ncomp+woffset), names(1:ncomp+woffset), types(1+ncomp+woffset:2*ncomp+woffset), state.matrix, ...
        pressures_matrix, sat_matrix, mob, rho_matrix, ...
        {}, comps_matrix, ...
        drivingForces);

    eqs(1+ncomp+woffset:2*ncomp+woffset) = eqsmat;
end

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
    [eqsW, namesW, typesW, state.wellSol, src] = model.insertWellEquations(eqs(1:ncomp+woffset), names(1:ncomp+woffset), ...
                                                     types(1:ncomp+woffset), wellSol0, wellSol, ...
                                                     wellvars, ...
                                                     wellMap, p, mob, rho, ...
                                                     {}, comps, ...
                                                     dt, opt); %#ok
end
% Wei
eqs = [eqsW(1:ncomp+woffset), eqs(1+ncomp+woffset:end),eqsW(1+ncomp+woffset:end)];
names = [namesW(1:ncomp+woffset), names(1+ncomp+woffset:end),namesW(1+ncomp+woffset:end)];
types = [typesW(1:ncomp+woffset), types(1+ncomp+woffset:end),typesW(1+ncomp+woffset:end)];

eq_offset = nwelleqs + 2*ncomp + 2*model.water;
for i = 1:ncomp
    eqs{i} = s.AccDiv(eqs{i}, fluxes{i}) + (TG{i}+TO{i}+Tdg{i});
    eqs{i+ncomp} = s.AccDiv(eqs{i+ncomp}, fluxes{i+ncomp}) - (TG{i}+TO{i}+Tdg{i});

    if isfield(model.rock, 'Mechanisms') && isfield(model.rock.Mechanisms,'diffusion')
        %Wei edit to add diffusion:
        Sg_poro = s.faceUpstr(upcg,sG.*model.rock.poro);
        Jg = s.Div(model.operators.T_diff{i}./model.rock.tau.*Sg_poro.*s.Grad(yM{i}.*rhoG));
        
        eqs{i} = eqs{i} - Jg;
    end
    
    ix = i + eq_offset;
    names{ix}= ['f_', mixture.names{i}];
    types{ix} = 'fugacity';
    eqs{ix} = (f_L{i}(twoPhase) - f_V{i}(twoPhase))/barsa;

    ix = i + eq_offset + ncomp;
    names{ix}= ['f_', mixture.names{i},'_matrix'];
    types{ix} = 'fugacity';
    eqs{ix} = (f_L_matrix{i}(twoPhase_matrix) - f_V_matrix{i}(twoPhase_matrix))/barsa;

    absent = state.components(twoPhase, i) <= 10*z_tol;
    if model.water
        absent = absent | pureWater(twoPhase);
    end
    if any(absent) && isa(eqs{ix}, 'ADI')
        eqs{ix}.val(absent) = 0;
    end
    
    absent_matrix = state.matrix.components(twoPhase_matrix, i) <= 10*z_tol;
    if model.water
        absent_matrix = absent_matrix | pureWater_matrix(twoPhase_matrix);
    end
    if any(absent_matrix) && isa(eqs{ix}, 'ADI')
        eqs{ix}.val(absent_matrix) = 0;
    end

    if opt.reduceToPressure
        C{ix - nwelleqs} = eqs{ix};
    end
end

if model.water
    eqs{ncomp+1} = s.AccDiv(eqs{ncomp+1}, rWvW);
end

cloix = eq_offset + 2*ncomp + 1;

if any(multiPhase)
    eqs{cloix} = sW(multiPhase) + sO(multiPhase) + sG(multiPhase) - 1;
else
    eqs{cloix} = [];
end
% Wei
if any(multiPhase_matrix)
    eqs{cloix+1} = sW_matrix(multiPhase_matrix) + sO_matrix(multiPhase_matrix) + sG_matrix(multiPhase_matrix) - 1;
else
    eqs{cloix+1} = [];
end

types{cloix} = 'saturation';
names{cloix} = 'volclosure';
% Wei
types{cloix+1} = 'saturation';
names{cloix+1} = 'volclosure_matrix';

if model.water
    eqs{ncomp+1} = eqs{ncomp+1}.*model.fluid.rhoWS;
end
if opt.reduceToPressure
    C{cloix - nwelleqs} = eqs{cloix};
    
    if model.water
        C{ncomp+1} = C{ncomp+1}.*model.fluid.rhoWS;
    end

    problem = PressureReducedLinearSystem(eqs, types, names, primaryVars, state, dt);
    problem.accumulationTerms = C;
    problem.model = model;
    problem.wellVarIndices = nc*(ncomp+model.water) + (1:nwellvar);
    problem.nwellvar = nwelleqs;
    problem.wellvars = wellvars;
    problem.wellvarNames = wellVarNames;
else
    if model.reduceLinearSystem
        problem = ReducedLinearizedSystem(eqs, types, names, primaryVars, state, dt);
    else
        problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
    end
end

if isa(problem, 'ReducedLinearizedSystem')
    % problem.keepNum = model.G.cells.num*(ncomp+model.water);
    problem.keepNum = nc*(2*ncomp+model.water) + nwellvar;
    problem.reorder = reorder;
end

problem.iterationNo = opt.iteration;
end


function [sO, sG] = setMinimums(model, state, sW, sO, sG, pureVapor, pureLiquid)
    stol = 1e-8;
    if model.water
        sT = sum(state.s, 2);
        if any(pureVapor)
            sG(pureVapor) = sT(pureVapor) - sW(pureVapor);
            if isa(sG, 'ADI')
                sG.val(pureVapor) = max(sG.val(pureVapor), stol);
            else
                sG(pureVapor) = max(sG(pureVapor), stol);
            end
        end

        if any(pureLiquid)
            sO(pureLiquid) = sT(pureLiquid) - sW(pureLiquid);
            if isa(sO, 'ADI')
                sO.val(pureLiquid) = max(sO.val(pureLiquid), stol);
            else
                sO(pureLiquid) = max(sO(pureLiquid), stol);
            end
        end
    end

end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

function [problem, state] = equationsCompositionalDPDP(state0, state, model, dt, drivingForces, varargin)
% Overall composition fully-implicit equations

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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
s_matrix = model.operators_matrix;
W = drivingForces.W;

fluid = model.fluid;
fluid_matrix = model.fluid_matrix;

% Properties at current timestep
[p, sW, z, temp, wellSol] = model.getProps(state, ...
                                                'pressure', 'sW', 'z', 'T', 'wellSol');
[p_matrix, sW_matrix, z_matrix] = model.getProps(state.matrix, ...
                                                'pressure', 'sW', 'z');

[p0, sW0, sO0, sG0, z0, temp0, wellSol0] = model.getProps(state0, ...
                                                         'pressure', 'sW', 'sO', 'sG', 'z', 'T', 'wellSol');
[p0_matrix, sW0_matrix, sO0_matrix, sG0_matrix, z0_matrix] = model.getProps(state0.matrix, ...
                                                         'pressure', 'sW', 'sO', 'sG', 'z');

z = ensureMinimumFraction(z);
z0 = ensureMinimumFraction(z0);
z = expandMatrixToCell(z);
z0 = expandMatrixToCell(z0);
z_matrix = ensureMinimumFraction(z_matrix);
z0_matrix = ensureMinimumFraction(z0_matrix);
z_matrix = expandMatrixToCell(z_matrix);
z0_matrix = expandMatrixToCell(z0_matrix);

[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

ncomp = numel(z);
cnames = model.EOSModel.getComponentNames();
cnames_matrix = cnames;
for i =1:ncomp
    cnames_matrix{i} = [cnames{i},'_matrix'];
end
nwellvar = sum(cellfun(@numel, wellVars));
if model.water
    [p, z{1:ncomp-1}, p_matrix, z_matrix{1:ncomp-1}, sW, wellVars{:}] = model.AutoDiffBackend.initVariablesAD(p, z{1:ncomp-1}, p_matrix, z_matrix{1:ncomp-1}, sW, wellVars{:});
    primaryVars = {'pressure', cnames{1:end-1},'pressure_matrix', cnames_matrix{1:end-1}, 'sW', wellVarNames{:}};
else
    [p, z{1:ncomp-1}, p_matrix, z_matrix{1:ncomp-1}, wellVars{:}] = model.AutoDiffBackend.initVariablesAD(p, z{1:ncomp-1}, p_matrix, z_matrix{1:ncomp-1}, wellVars{:});
    primaryVars = {'pressure', cnames{1:end-1},'pressure_matrix', cnames_matrix{1:end-1}, wellVarNames{:}};
end
% Property pressure different from flow potential
p_prop = opt.propsPressure;
otherPropPressure = ~isempty(p_prop);
if ~otherPropPressure
    p_prop = p;
end
p_prop_matrix = p_matrix;

z{end} = ones(numel(value(p)),1); z_matrix{end} = z{end} ;
for i = 1:(ncomp-1)
    z{end} = z{end} - z{i};
    z_matrix{end} = z_matrix{end} - z_matrix{i};
end

[xM,  yM,  sO,  sG,  rhoO,  rhoG, muO, muG] = model.computeTwoPhaseFlowProps(state, p_prop, temp, z);
[xM0, yM0, ~, ~, rhoO0, rhoG0] = model.computeTwoPhaseFlowProps(state0, p0, temp0, z0);
[xM_matrix,  yM_matrix,  sO_matrix,  sG_matrix,  rhoO_matrix,  rhoG_matrix, muO_matrix, muG_matrix] = model.computeTwoPhaseFlowProps(state.matrix, p_prop_matrix, temp, z_matrix);
[xM0_matrix, yM0_matrix, ~, ~, rhoO0_matrix, rhoG0_matrix] = model.computeTwoPhaseFlowProps(state0.matrix, p0_matrix, temp0, z0_matrix);

state.rho = [value(rhoO), value(rhoG)];
state.matrix.rho = [value(rhoO_matrix), value(rhoG_matrix)];
% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p_prop, p0);
[pvMult_matrix, transMult_matrix, mobMult_matrix, pvMult0_matrix] = getMultipliers(model.fluid_matrix, p_prop_matrix, p0_matrix);

if model.water
    sO = (1-sW).*sO;
    sG = (1-sW).*sG;
    sat = {sW, sO, sG};
    
    [krW, krO, krG] = model.evaluateRelPerm(sat);
    krW = mobMult.*krW;
    % Wei: matrix
    sO_matrix = (1-sW_matrix).*sO_matrix;
    sG_matrix = (1-sW_matrix).*sG_matrix;
    sat_matrix = {sW_matrix, sO_matrix, sG_matrix};
    
    [krW_matrix, krO_matrix, krG_matrix] = model.evaluateRelPerm(sat_matrix);
    krW_matrix = mobMult_matrix.*krW_matrix;
else
    sat = {sO, sG};
    sat_matrix = {sO_matrix, sG_matrix};

    if ~isfield(model.fluid, 'relPermScal')
        [krO, krG] = model.evaluateRelPerm(sat,'medium','fracture');
        [krO_matrix, krG_matrix] = model.evaluateRelPerm(sat_matrix,'medium','matrix');
    else
        % Wei edit Relative Permeability Scaling, but some function is missing
        kr = eCPARelativePermeabilityScaling(model, state);
        kr_matrix = eCPARelativePermeabilityScaling(model, state.matrix);
        [krO, krG] = deal(kr{:});
        [krO_matrix, krG_matrix] = deal(kr_matrix{:});
    end
end

krO = mobMult.*krO;
krG = mobMult.*krG;
% Wei: Modifiy relperm by mobility multiplier (if any)
krO_matrix = mobMult_matrix.*krO_matrix;
krG_matrix = mobMult_matrix.*krG_matrix;

% Compute transmissibility
T = s.T.*transMult;
T_matrix = s_matrix.T.*transMult_matrix;

% Gravity gradient per face
gdz = model.getGravityGradient();
gdz_matrix = model.getGravityGradient();

% Oil flux
rhoOf  = s.faceAvg(sO.*rhoO)./max(s.faceAvg(sO), 1e-8);
mobO   = krO./muO;
dpO    = s.Grad(p) - rhoOf.*gdz;
upco  = (value(dpO)<=0);
vO = -s.faceUpstr(upco, mobO).*T.*dpO;

rhoOf_matrix  = s_matrix.faceAvg(sO_matrix.*rhoO_matrix)./max(s_matrix.faceAvg(sO_matrix), 1e-8);
mobO_matrix   = krO_matrix./muO_matrix;
dpO_matrix    = s_matrix.Grad(p_matrix) - rhoOf_matrix.*gdz_matrix;
upco_matrix  = (value(dpO_matrix)<=0);
vO_matrix = -s_matrix.faceUpstr(upco_matrix, mobO_matrix).*T_matrix.*dpO_matrix;

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
% Wei: matrix
rhoGf_matrix  = s_matrix.faceAvg(sG_matrix.*rhoG_matrix)./max(s_matrix.faceAvg(sG_matrix), 1e-8);
mobG_matrix   = krG_matrix./muG_matrix;
if isfield(fluid_matrix, 'pcOG')
    pcOG_matrix  = fluid_matrix.pcOG(sG);
    pG_matrix = p_matrix + pcOG_matrix;
else
    pG_matrix = p_matrix;
end
dpG_matrix  = s_matrix.Grad(pG_matrix) - rhoGf_matrix.*gdz_matrix;
upcg_matrix = (value(dpG_matrix)<=0);
vG_matrix   = -s_matrix.faceUpstr(upcg_matrix, mobG_matrix).*T_matrix.*dpG_matrix;

% mass flux
rOvO = s.faceUpstr(upco, rhoO).*vO;
rGvG = s.faceUpstr(upcg, rhoG).*vG;
rOvO_matrix = s_matrix.faceUpstr(upco_matrix, rhoO_matrix).*vO_matrix;
rGvG_matrix = s_matrix.faceUpstr(upcg_matrix, rhoG_matrix).*vG_matrix;

bO = rhoO./fluid.rhoOS;
bG = rhoG./fluid.rhoGS;
bO_matrix = rhoO_matrix./fluid_matrix.rhoOS;
bG_matrix = rhoG_matrix./fluid_matrix.rhoGS;

% EQUATIONS -----------------------------------------------------------
% Wei: water equation + 2n component equations
[eqs, fluxes, types, names] = deal(cell(1, 2*ncomp + model.water));

woffset = model.water;
if model.water
    % Water flux
    muW = fluid.muW(p_prop);
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

    % Wei: Water flux of matrix
    muW_matrix = fluid_matrix.muW(p_prop_matrix);
    bW_matrix  = fluid_matrix.bW(p_prop_matrix);
    rhoW_matrix= bW_matrix.*fluid_matrix.rhoWS;
    bW0_matrix = fluid_matrix.bW(p0_matrix);

    rhoWf_matrix  = s_matrix.faceAvg(rhoW_matrix);
    mobW_matrix   = krW_matrix./muW_matrix;
    if isfield(fluid_matrix, 'pcOW')
        pcOW_matrix  = fluid_matrix.pcOW(sW_matrix);
        pW_matrix = p_matrix - pcOW_matrix;
    else
        pW_matrix = p_matrix;
    end
    dpW_matrix    = s_matrix.Grad(pW_matrix) - rhoWf_matrix.*gdz_matrix;
    upcw_matrix  = (value(dpW_matrix)<=0);
    vW_matrix = -s_matrix.faceUpstr(upcw_matrix, mobW_matrix).*T_matrix.*dpW_matrix;
    rWvW_matrix = s_matrix.faceUpstr(upcw_matrix, bW_matrix).*vW_matrix;
    
    eqs{ncomp+1+woffset} = (s_matrix.pv/dt).*( bW_matrix.*pvMult_matrix.*sW_matrix - bW0_matrix.*pvMult0_matrix.*sW0_matrix);
    names{ncomp+1+woffset} = 'water_matrix';
    types{ncomp+1+woffset} = 'cell';
    fluxes{ncomp+1+woffset} = rWvW_matrix;
else
    [vW, mobW, upcw, bW, rhoW] = deal([]);
    [vW_matrix, mobW_matrix, upcw_matrix, bW_matrix, rhoW_matrix] = deal([]);
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
krGt{i} = krG.*yM{i}.*rhoG./muG.*dpsiG + krG_matrix.*yM_matrix{i}.*rhoG_matrix./muG_matrix.*(~dpsiG);
krOt{i} = krO.*xM{i}.*rhoO./muO.*dpsiO + krO_matrix.*xM_matrix{i}.*rhoO_matrix./muO_matrix.*(~dpsiO);

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

%
if model.outputFluxes
    % Wei: save internal fluxes
    state.matrix = model.storeFluxes(state.matrix, vW_matrix, vO_matrix, vG_matrix);
    state = model.storeFluxes(state, vW, vO, vG);
end
if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, bG);
    state = model.storeMobilities(state, mobW, mobO, mobG);
    state = model.storeUpstreamIndices(state, upcw, upco, upcg);
    state = model.storeDensities(state, rhoW, rhoO, rhoG);

    state.matrix = model.storebfactors(state.matrix, bW_matrix, bO_matrix, bG_matrix);
    state.matrix = model.storeMobilities(state.matrix, mobW_matrix, mobO_matrix, mobG_matrix);
    state.matrix = model.storeUpstreamIndices(state.matrix, upcw_matrix, upco_matrix, upcg_matrix);
    state.matrix = model.storeDensities(state.matrix, rhoW_matrix, rhoO_matrix, rhoG_matrix);
end

if isfield(model.rock, 'Mechanisms') && isfield(model.rock.Mechanisms,'sorption')
    sumIso = yM_matrix{1}./model.rock.isothermP(1);
    sumIso0 = yM0_matrix{1}./model.rock.isothermP(1);
    for ii=2:numel(yM_matrix)
        sumIso = sumIso + yM_matrix{ii}./model.rock.isothermP(ii);
        sumIso0 = sumIso0 + yM0_matrix{ii}./model.rock.isothermP(ii);
    end
end

for i = 1:ncomp
    names{i+woffset} = cnames{i};
    types{i+woffset} = 'cell';

    eqs{i+woffset} = (s.pv/dt).*( ...
                    rhoO.*pvMult.*sO.*xM{i} - rhoO0.*pvMult0.*sO0.*xM0{i} + ...
                    rhoG.*pvMult.*sG.*yM{i} - rhoG0.*pvMult0.*sG0.*yM0{i});
    fluxes{i+woffset} = rOvO.*s.faceUpstr(upco, xM{i}) + rGvG.*s.faceUpstr(upcg, yM{i});
    % Wei: matrix
    names{i+woffset+ncomp} = cnames_matrix{i};
    types{i+woffset+ncomp} = 'cell';

    if isfield(model.rock, 'Mechanisms') && isfield(model.rock.Mechanisms,'sorption')
        eqs{i+woffset+ncomp} = (1/dt).*( ...
            s_matrix.pv.*rhoO_matrix.*pvMult_matrix.*sO_matrix.*xM_matrix{i} - s_matrix.pv.*rhoO0_matrix.*pvMult0_matrix.*sO0_matrix.*xM0_matrix{i} + ...
            s_matrix.pv.*rhoG_matrix.*pvMult_matrix.*sG_matrix.*yM_matrix{i} - s_matrix.pv.*rhoG0_matrix.*pvMult0_matrix.*sG0_matrix.*yM0_matrix{i} +...
            (s_matrix.gv.*model.rock.isothermRho(i)./model.rock.isothermP(i)).*...
            ((yM_matrix{i}.*p_matrix)./(1+p_matrix.*sumIso) - (yM0_matrix{i}.*p0_matrix)./(1+p0_matrix.*sumIso0)) );
    else
        eqs{i+woffset+ncomp} = (s_matrix.pv/dt).*( ...
            rhoO_matrix.*pvMult_matrix.*sO_matrix.*xM_matrix{i} - rhoO0_matrix.*pvMult0_matrix.*sO0_matrix.*xM0_matrix{i} + ...
            rhoG_matrix.*pvMult_matrix.*sG_matrix.*yM_matrix{i} - rhoG0_matrix.*pvMult0_matrix.*sG0_matrix.*yM0_matrix{i});
    end
    fluxes{i+woffset+ncomp} = rOvO_matrix.*s_matrix.faceUpstr(upco_matrix, xM_matrix{i}) + rGvG_matrix.*s_matrix.faceUpstr(upcg_matrix, yM_matrix{i});
end

% Finally, add in and setup well equations
if model.water
    rho = {rhoW, rhoO, rhoG};
    mob = {mobW, mobO, mobG};
    pressures = {pW, p, pG};

    %Wei
    rho_matrix = {rhoW_matrix, rhoO_matrix, rhoG_matrix};
    mob_matrix = {mobW_matrix, mobO_matrix, mobG_matrix};
    pressures_matrix = {pW_matrix, p_matrix, pG_matrix};
else
    rho = {rhoO, rhoG};
    mob = {mobO, mobG};
    pressures = {p, pG};

    %Wei
    rho_matrix = {rhoO_matrix, rhoG_matrix};
    mob_matrix = {mobO_matrix, mobG_matrix};
    pressures_matrix = {p_matrix, pG_matrix};
end
comps = cellfun(@(x, y) {x, y}, xM, yM, 'UniformOutput', false);
comps_matrix = cellfun(@(x, y) {x, y}, xM_matrix, yM_matrix, 'UniformOutput', false);

[eqsBC, state] = model.addBoundaryConditionsAndSources(eqs(1:ncomp+woffset), names(1:ncomp+woffset), types(1:ncomp+woffset), state, ...
                                                     pressures, sat, mob, rho, ...
                                                     {}, comps, ...
                                                     drivingForces);
% Wei
eqs(1:ncomp+woffset) = eqsBC;

if ~isempty(drivingForces.bc)
    % Compute the contribution from BC while avoiding the capillary jump

    % Store the fracture operators and BC values
    frop = model.operators;
    frbc = drivingForces.bc.value;

    % Set up the boundary conditions using the matrix operators and BC values
    model.operators = model.operators_matrix;
    drivingForces.bc.value = drivingForces.bc.value_matrix;

    [eqsmat, state.matrix] = addBoundaryConditionsAndSources(model, ...
        eqs(1+ncomp+woffset:2*ncomp+woffset), names(1:ncomp+woffset), types(1+ncomp+woffset:2*ncomp+woffset), state.matrix, ...
        pressures_matrix, sat_matrix, mob_matrix, rho_matrix, ...
        {}, comps_matrix, ...
        drivingForces);

    eqs(1+ncomp+woffset:2*ncomp+woffset) = eqsmat;

    % Restore the fracture operators and BC values
    model.operators = frop;
    drivingForces.bc.value = frbc;
end

if opt.staticWells
    compSrc = vertcat(wellSol.components);
    for i = 1:ncomp
        wc = vertcat(W.cells);
        eqs{i+woffset}(wc) = eqs{i+woffset}(wc) - compSrc(:, i);
    end
else
    [eqsW, namesW, typesW, state.wellSol] = model.insertWellEquations(eqs(1:ncomp+woffset), names(1:ncomp+woffset), ...
                                                     types(1:ncomp+woffset), wellSol0, wellSol, ...
                                                     wellVars, ...
                                                     wellMap, p, mob, rho, ...
                                                     {}, comps, ...
                                                     dt, opt);
end
% Wei
eqs = [eqsW(1:ncomp+woffset), eqs(1+ncomp+woffset:end),eqsW(1+ncomp+woffset:end)];
names = [namesW(1:ncomp+woffset), names(1+ncomp+woffset:end),namesW(1+ncomp+woffset:end)];
types = [typesW(1:ncomp+woffset), types(1+ncomp+woffset:end),typesW(1+ncomp+woffset:end)];

conservationindices = 1:(2*ncomp+woffset);
acc = eqs(conservationindices);
for i = 1:(ncomp+woffset)
    eqs{i} = s.AccDiv(eqs{i}, fluxes{i}) + (TG{i}+TO{i}+Tdg{i});
    % Wei
    eqs{i+ncomp+woffset} = s.AccDiv(eqs{i+ncomp+woffset}, fluxes{i+ncomp+woffset}) - (TG{i}+TO{i}+Tdg{i});
    if isfield(model.rock, 'Mechanisms') && isfield(model.rock.Mechanisms,'diffusion')
        %Wei edit to add diffusion:
        Sg_poro = s.faceUpstr(upcg,sG.*model.rock.poro);
        Jg = s.Div(model.operators.T_diff{i}./model.rock.tau.*Sg_poro.*s.Grad(yM{i}.*rhoG));
        
        Sg_poro_matrix = s_matrix.faceUpstr(upcg_matrix,sG_matrix.*model.rock_matrix.poro);
        Jg_matrix = s_matrix.Div(model.operators.T_diff{i}./model.rock.tau.*Sg_poro_matrix.*s_matrix.Grad(yM_matrix{i}.*rhoG_matrix));

        eqs{i} = eqs{i} - Jg;
        eqs{i+ncomp} = eqs{i+ncomp} - Jg_matrix;
    end

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
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

function [problem, state] = pressureEquationBlackOilPolymer(state0, state, model, dt, drivingForces, varargin)

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'staticWells',  false, ...
             'propsPressure', [], ...
             'iteration', -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.W;
op = model.operators;
fluid = model.fluid;

disgas = model.disgas;
vapoil = model.vapoil;

% Polymer shear thinning/thickening
usingShear = isfield(fluid, 'plyshearMult');

% Properties at current timestep
[p, sW, sG, c, cmax, rs, rv, wellSol] = model.getProps(state, ...
                                'pressure', 'water', 'gas', 'polymer', ...
                                'polymermax', 'rs', 'rv', 'wellSol');
% Properties at previous timestep
[p0, sW0, sG0, c0, cmax0, rs0, rv0, wellSol0] = model.getProps(state0, ...
                                'pressure', 'water', 'gas', 'polymer', ...
                                'polymermax', 'rs', 'rv', 'wellSol');


[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

%Initialization of independent variables ----------------------------------
st  = model.getCellStatusVO(state,  1-sW-sG,   sW,  sG);
st0 = model.getCellStatusVO(state0, 1-sW0-sG0, sW0, sG0);
p_prop = opt.propsPressure;
otherPropPressure = ~isempty(p_prop);
if ~opt.resOnly,
    if ~opt.reverseMode,
        % define primary varible x and initialize
        x = st{1}.*rs + st{2}.*rv + st{3}.*sG;

        [p, wellVars{:}] = model.AutoDiffBackend.initVariablesAD(p, wellVars{:});
        if ~otherPropPressure
            p_prop = p;
        end
        % define sG, rs and rv in terms of x
        sG = st{2}.*(1-sW) + st{3}.*x;
        if disgas
            rsSat = fluid.rsSat(p_prop);
            rs = (~st{1}).*rsSat + st{1}.*x;
        else % otherwise rs = rsSat = const
            rsSat = rs;
        end
        if vapoil
            rvSat = fluid.rvSat(p_prop);
            rv = (~st{2}).*rvSat + st{2}.*x;
        else % otherwise rv = rvSat = const
            rvSat = rv;
        end
    else
        assert(0, 'Backwards solver not supported for splitting');
    end
else % resOnly-case compute rsSat and rvSat for use in well eqs
    if isempty(p_prop)
        p_prop = p;
    end
    if disgas, rsSat = fluid.rsSat(p_prop); else rsSat = rs; end
    if vapoil, rvSat = fluid.rvSat(p_prop); else rvSat = rv; end
end
sO  = 1- sW  - sG;
sO0 = 1- sW0 - sG0;

primaryVars = {'pressure', wellVarNames{:}};

% FLIUD PROPERTIES ---------------------------------------------------
[krW, krO, krG] = model.evaluateRelPerm({sW, sO, sG});

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p_prop, p0);

% Modifiy relperm by mobility multiplier (if any)
krW = mobMult.*krW; krO = mobMult.*krO; krG = mobMult.*krG;

% Compute transmissibility
T = op.T.*transMult;

% Gravity gradient per face
gdz = model.getGravityGradient();

% Evaluate adsorptions
ads  = effads(c, cmax, model);
ads0  = effads(c0, cmax0, model);

% Evaluate water and polymer properties
[vW, vP, bW, muWMult, mobW, mobP, rhoW, pW, upcw, a, dpW] = ...
    getFluxAndPropsWaterPolymer_BO(model, p_prop, sW, c, ads, ...
    krW, T, gdz); %, 'shear', ~otherPropPressure);
bW0 = fluid.bW(p0);

% Evaluate oil properties
[vO, bO, mobO, rhoO, pO, upco, dpO] = getFluxAndPropsOil_BO(model, p_prop, sO, krO, T, gdz, rs, ~st{1});
bO0 = getbO_BO(model, p0, rs0, ~st0{1});

% Evaluate gas properties
[vG, bG, mobG, rhoG, pG, upcg, dpG] = getFluxAndPropsGas_BO(model, p_prop, sG, krG, T, gdz, rv, ~st{2});
bG0 = getbG_BO(model, p0, rv0, ~st0{2});

if otherPropPressure
    % We have used a different pressure for property evaluation, undo the
    % effects of this on the fluxes.
    dp_diff = op.Grad(p) - op.Grad(p_prop);
    
    vW = -op.faceUpstr(upcw, mobW).*op.T.*(dpW + dp_diff);
    vO = -op.faceUpstr(upco, mobO).*op.T.*(dpO + dp_diff);
    vG = -op.faceUpstr(upcg, mobG).*op.T.*(dpG + dp_diff);
    vP = -op.faceUpstr(upcw, mobP).*op.T.*(dpW + dp_diff);
    
        % If other property pressure is used, then the shear thinning is
    % computed here, and not above
    if usingShear
        poroFace  = op.faceAvg(model.rock.poro);
        faceArea  = model.G.faces.areas(op.internalConn);
        Vw        = vW./(poroFace .* faceArea); % water velocity
        muWMultf  = op.faceUpstr(upcw, muWMult);
        [shearMult, ~, shearReport] = ...
            getPolymerShearMultiplier(model, Vw, muWMultf);
        vW        = vW .* shearMult;
        vP        = vP .* shearMult;
        extraOutput.shearMult   = shearMult;
        extraOutput.shearReport = shearReport;
    end
end

% These are needed in transport solver, so we output them regardless of
% any flags set in the model.
state = model.storeFluxes(state, vW, vO, vG);
state = model.storeUpstreamIndices(state, upcw, upco, upcg);
if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, bG);
    state = model.storeMobilities(state, mobW, mobO, mobG, mobP);
end

if usingShear
    % The shear multipliers are stored for use in the transport solver
    state = model.storeShearMultiplier(state, extraOutput.shearMult);
end

% EQUATIONS -----------------------------------------------------------

% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.
bOvO = op.faceUpstr(upco, bO).*vO;
bWvW = op.faceUpstr(upcw, bW).*vW;
bGvG = op.faceUpstr(upcg, bG).*vG;
bWvP = op.faceUpstr(upcw, bW).*vP;

% The first equation is the conservation of the water phase. This equation is
% straightforward, as water is assumed to remain in the aqua phase in the
% black oil model.
water = (op.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + op.Div(bWvW);

% Same goes for the polymer equation
poro = model.rock.poro;
polymer = (op.pv.*(1-fluid.dps)/dt).*(pvMult.*bW.*sW.*c - ...
   pvMult0.*bW0.*sW0.*c0) + (op.pv/dt).* ...
   (fluid.rhoR.*((1-poro)./poro).*(ads-ads0) ) + op.Div(bWvP);

if ~opt.resOnly
    epsilon = 1.e-8;
    % the first way is based on the diagonal values of the resulting
    % Jacobian matrix
    if isa(polymer, 'ADI')     
        eps = sqrt(epsilon)*mean(abs(diag(polymer.jac{3})));
        % bad marks the cells prolematic in evaluating Jacobian
        bad = abs(diag(polymer.jac{3})) < eps;
        % the other way is to choose based on the water saturation
        polymer(bad) = c(bad);
    end
else
    assert(0, 'Backwards solver not supported for splitting');
end

% Second equation: mass conservation equation for the oil phase at surface
% conditions. This is any liquid oil at reservoir conditions, as well as
% any oil dissolved into the gas phase (if the model has vapoil enabled).
if model.vapoil
    % The model allows oil to vaporize into the gas phase. The conservation
    % equation for oil must then include the fraction present in the gas
    % phase.
    rvbGvG = op.faceUpstr(upcg, rv).*bGvG;
    % Final equation
    oil = (op.pv/dt).*( pvMult.* (bO.* sO  + rv.* bG.* sG) - ...
        pvMult0.*(bO0.*sO0 + rv0.*bG0.*sG0) ) + ...
        op.Div(bOvO + rvbGvG);
else
    oil = (op.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) + op.Div(bOvO);
end

% Conservation of mass for gas. Again, we have two cases depending on
% whether the model allows us to dissolve the gas phase into the oil phase.
if model.disgas
    % The gas transported in the oil phase.
    rsbOvO = op.faceUpstr(upco, rs).*bOvO;
    
    gas = (op.pv/dt).*( pvMult.* (bG.* sG  + rs.* bO.* sO) - ...
        pvMult0.*(bG0.*sG0 + rs0.*bO0.*sO0 ) ) + ...
        op.Div(bGvG + rsbOvO);
else
    gas = (op.pv/dt).*( pvMult.*bG.*sG - pvMult0.*bG0.*sG0 ) + op.Div(bGvG);
end
eqs = {water, oil, gas, polymer};
names = {'water', 'oil', 'gas', 'polymer'};
types = {'cell', 'cell', 'cell', 'cell'};

rho = {rhoW, rhoO, rhoG};
mob = {mobW, mobO, mobG};
sat = {sW, sO, sG};
dissolved = model.getDissolutionMatrix(rs, rv);

[eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                             {pW, p, pG}, sat, mob, rho, ...
                                             dissolved, {c}, ...
                                             drivingForces);

% Finally, add in and setup well equations
[eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, ...
    types, wellSol0, wellSol, wellVars, wellMap, p, mob, rho, dissolved, {c}, dt, opt);

% Create actual pressure equation
cfac = 1./(1 - disgas*vapoil*rs.*rv);

a_w = 1./bW;
a_o = cfac.*(1./bO - disgas*rs./bG);
a_g = cfac.*(1./bG - vapoil*rv./bO);

water = eqs{1};
oil   = eqs{2};
gas   = eqs{3};

eqs{1} = (dt./op.pv).*(oil.*a_o + water.*a_w + gas.*a_g);
names{1} = 'pressure';
types{1} = 'cell';

% Strip phase equations
eqs = eqs([1, 5:end]);
names = names([1, 5:end]);
types = types([1, 5:end]);

state.timestep = dt;
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end

%--------------------------------------------------------------------------
% Helper functions
%--------------------------------------------------------------------------

% Effective adsorption, depending of desorption or not
function y = effads(c, cmax, model)
   if model.fluid.adsInx == 2
      y = model.fluid.ads(max(c, cmax));
   else
      y = model.fluid.ads(c);
   end
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


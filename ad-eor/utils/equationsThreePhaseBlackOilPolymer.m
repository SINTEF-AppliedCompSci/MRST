function [problem, state] = equationsThreePhaseBlackOilPolymer(state0, state, model, dt, drivingForces, varargin)

%
%
% SYNOPSIS:
%   function [problem, state] = equationsThreePhaseBlackOilPolymer(state0, state, model, dt, drivingForces, varargin)
%
% DESCRIPTION: 
%   Assemble the linearized equations for a blackoil system,
%   computing both the residuals and the Jacobians. Returns the result as an
%   instance of the class LinearizedProblem which can be solved using instances of
%   LinearSolverAD.
%
%   A description of the modeling equations can be found in the directory
%   ad-eor/docs.
%
% PARAMETERS:
%   state0        - State at previous times-step
%   state         - State at current time-step
%   model         - Model instance
%   dt            - time-step
%   drivingForces - Driving forces (boundary conditions, wells, ...)
%   varargin      - optional parameters
%
% RETURNS:
%   problem - Instance of LinearizedProblem
%   state   - Updated state variable (fluxes, mobilities and more can be
%             stored, the wellSol structure is also updated in case of control switching)
%
% EXAMPLE:
%
% SEE ALSO: LinearizedProblem, LinearSolverAD, equationsOilWater, OilWaterPolymerModel
%

% Get linearized problem for oil/water/polymer system with black oil-style
% properties

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

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'iteration', -1);

opt = merge_options(opt, varargin{:});

% Shorter names for some commonly used parts of the model and forces.
s = model.operators;
f = model.fluid;
G = model.G;
W = drivingForces.W;
% Currently we do not support senario without wells.
assert(isempty(drivingForces.bc) && isempty(drivingForces.src));

% Properties at current timestep
[p, sW, sG, rs, rv, cp, cpmax, wellSol] = model.getProps(state, ...
   'pressure', 'water', 'gas', 'rs', 'rv', 'polymer', 'polymermax', 'wellsol');

% Properties at previous timestep
[p0, sW0, sG0, rs0, rv0, cp0, cpmax0, wellSol0] = model.getProps(state0, ...
   'pressure', 'water', 'gas', 'rs', 'rv', 'polymer', 'polymermax', 'wellsol');

[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

%Initialization of primary variables ----------------------------------
st  = model.getCellStatusVO(state,  1-sW-sG,   sW,  sG);
st0 = model.getCellStatusVO(state0, 1-sW0-sG0, sW0, sG0);

if model.disgas || model.vapoil
    % X is either Rs, Rv or Sg, depending on each cell's saturation status
    x = st{1}.*rs + st{2}.*rv + st{3}.*sG;
    gvar = 'x';
else
    x = sG;
    gvar = 'sG';
end

if ~opt.resOnly
    if ~opt.reverseMode
        % define primary varible x and initialize
        [p, sW, x, cp, wellVars{:}] = ...
            initVariablesADI(p, sW, x, cp, wellVars{:});
    else
        x0 = st0{1}.*rs0 + st0{2}.*rv0 + st0{3}.*sG0;
        % Set initial gradient to zero
        wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
        [p0, sW0, x0, cp0, wellVars0{:}] = ...
            initVariablesADI(p0, sW0, x0, cp0, wellVars0{:}); %#ok
        clear zw;
        [sG0, rs0, rv0] = calculateHydrocarbonsFromStatusBO(model, st0, 1-sW, x0, rs0, rv0, p0);
    end
end

if ~opt.reverseMode
    % Compute values from status flags. If we are in reverse mode, these
    % values have already converged in the forward simulation.
    [sG, rs, rv, rsSat, rvSat] = calculateHydrocarbonsFromStatusBO(model, st, 1-sW, x, rs, rv, p);
end

% We will solve for pressure, water and gas saturation (oil saturation
% follows via the definition of saturations), polymer concentration and well rates + bhp.
primaryVars = {'pressure', 'sW', gvar, 'polymer', wellVarNames{:}};

% Evaluate relative permeability
sO  = 1 - sW  - sG;
sO0 = 1 - sW0 - sG0;
[krW, krO, krG] = model.evaluateRelPerm({sW, sO, sG});

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

% Modifiy relperm by mobility multiplier (if any)
krW = mobMult.*krW; krO = mobMult.*krO; krG = mobMult.*krG;

% Compute transmissibility
T = s.T.*transMult;

% Gravity gradient per face
gdz = model.getGravityGradient();

% Evaluate water and polymer properties
ads  = effads(cp, cpmax, f);
ads0 = effads(cp0, cpmax0, f);
[vW, vP, bW, muWMult, mobW, mobP, rhoW, pW, upcw, a] = ...
    getFluxAndPropsWaterPolymer_BO(model, p, sW, cp, ads, ...
    krW, T, gdz);
bW0 = f.bW(p0);


% Evaluate oil properties
[vO, bO, mobO, rhoO, p, upco] = getFluxAndPropsOil_BO(model, p, sO, krO, T, gdz, rs, ~st{1});
bO0 = getbO_BO(model, p0, rs0, ~st0{1});

% Evaluate gas properties
bG0 = getbG_BO(model, p0, rv0, ~st0{2});
[vG, bG, mobG, rhoG, pG, upcg] = getFluxAndPropsGas_BO(model, p, sG, krG, T, gdz, rv, ~st{2});

if model.usingShear || model.usingShearLog || model.usingShearLogshrate
    % calculate well perforation rates :
    if ~isempty(W)
        if ~opt.reverseMode
            wc    = vertcat(W.cells);
            perf2well = getPerforationToWellMapping(W);
            sgn = vertcat(W.sign);
            
            wc_inj = wc(sgn(perf2well) > 0);
            cpw = cp(wc_inj);
            
            mob = {mobW, mobO, mobG};
            muWMultW = muWMult(wc_inj);
            muWFullyMixed = model.fluid.muWMult(cpw);

            mob{1}(wc_inj) = mob{1}(wc_inj) ./ muWFullyMixed .* muWMultW;
            rho = {rhoW, rhoO, rhoG};
            
            dissolved = model.getDissolutionMatrix(rs, rv);

            
            [src, wellsys, state.wellSol] = ...
            model.FacilityModel.getWellContributions(wellSol0, wellSol, wellVars, ...
                                    wellMap, p, mob, rho, dissolved, {cp}, ...
                                    dt, opt.iteration);

        else
            error('not supported yet!');
        end
    else
        error('The polymer model does not support scenarios without wells now!');
    end

    % s = model.operators;  % The previous s was overwritten with saturations.
    poro =  s.pv./G.cells.volumes;
    poroFace = s.faceAvg(poro);
    faceA = G.faces.areas(s.internalConn);

    % Bw * Fw should be flux
    Vw = vW./(poroFace .* faceA);

    % Using the upstreamed viscosity multiplier due to PLYVISC
    muWMultf = s.faceUpstr(upcw, muWMult);

    wc = vertcat(W.cells);
    muWMultW = muWMult(wc);

    % We assume the viscosity multiplier should be consistent with current
    % way in handling the injection mobility, while the assumption is not
    % verfied with any tests yet due to lack of the reference result.

    [~, wciPoly, iInxW] = getWellPolymer(W);
    cpw = cp(wc);
    muWMultW(iInxW) = model.fluid.muWMult(cpw(iInxW));

    % Maybe should also apply this for PRODUCTION wells.
    muWMultW((iInxW(wciPoly==0))) = 1;


    % The water flux for the wells.
    cqs = vertcat(state.wellSol.cqs);
    fluxWaterWell = value(cqs(:, 1));

    poroW = poro(wc);

    % the thickness of the well perforations in the cell
    welldir = { W.dir };
    i = cellfun('prodofsize', welldir) == 1;
    welldir(i) = arrayfun(@(w) repmat(w.dir, [ numel(w.cells), 1 ]), ...
                          W(i), 'UniformOutput', false);
    welldir = vertcat(welldir{:});
    [dx, dy, dz] = cellDims(G, wc);
    thicknessWell = dz;
    thicknessWell(welldir == 'Y') = dy(welldir == 'Y');
    thicknessWell(welldir == 'X') = dx(welldir == 'X');

    % For the wells
    % The water velocity is computed at the reprensentative radius rR.
    if ~isfield(W, 'rR')
        error('The representative radius of the well is not initialized');
    end
    rR = vertcat(W.rR);

    VwW = bW(wc).*fluxWaterWell./(poroW .* rR .* thicknessWell * 2 * pi);

    muWMultW = value(muWMultW);
    VwW = value(VwW);
    muWMultf = value(muWMultf);
    Vw = value(Vw);


    if model.usingShearLogshrate
        % calculating the shear rate based on the velocity
        if ~opt.resOnly
            krwF = s.faceUpstr(upcw, krW.val);
            swF = s.faceUpstr(upcw, sW.val);
        end

        if opt.resOnly
            krwF = s.faceUpstr(upcw, krW);
            swF = s.faceUpstr(upcw, sW);
        end

        permF = s.T ./faceA;
        temp = permF.*swF.*krwF;

        index = find(abs(Vw) > 0.);
        Vw(index) = 4.8 * Vw(index).*sqrt(poroFace(index) ./temp(index));

        % calculating the shear rate for the wells
        rW = vertcat(W.r);
        VwW = 4.8 * VwW ./(2*rW);
    end

    if model.usingShear
        shearMultf = computeShearMult(model.fluid, abs(Vw), muWMultf);
        shearMultW = computeShearMult(model.fluid, abs(VwW), muWMultW);
    end

    if model.usingShearLog || model.usingShearLogshrate
        shearMultf = computeShearMultLog(model.fluid, abs(Vw), muWMultf);
        shearMultW = computeShearMultLog(model.fluid, abs(VwW), muWMultW);
    end

    vW = vW ./ shearMultf;
    vP = vP ./ shearMultf;
end

% EQUATIONS -----------------------------------------------------------

% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.
bOvO = s.faceUpstr(upco, bO).*vO;
bWvW = s.faceUpstr(upcw, bW).*vW;
bGvG = s.faceUpstr(upcg, bG).*vG;
bWvP = s.faceUpstr(upcw, bW).*vP;

% Store fluxes / properties for debugging / plotting, if requested.
if model.outputFluxes
    state = model.storeFluxes(state, vW, vO, vG);
end
if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, bG);
    state = model.storeMobilities(state, mobW, mobO, mobG);
    state = model.storeUpstreamIndices(state, upcw, upco, upcg);
end
% EQUATIONS ---------------------------------------------------------------

% The first equation is the conservation of the water phase. This equation is
% straightforward, as water is assumed to remain in the aqua phase in the
% black oil model.
water = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + s.Div(bWvW);

% Second equation: mass conservation equation for the oil phase at surface
% conditions. This is any liquid oil at reservoir conditions, as well as
% any oil dissolved into the gas phase (if the model has vapoil enabled).
if model.vapoil
    % The model allows oil to vaporize into the gas phase. The conservation
    % equation for oil must then include the fraction present in the gas
    % phase.
    rvbGvG = s.faceUpstr(upcg, rv).*bGvG;
    % Final equation
    oil = (s.pv/dt).*( pvMult.* (bO.* sO  + rv.* bG.* sG) - ...
        pvMult0.*(bO0.*sO0 + rv0.*bG0.*sG0) ) + ...
        s.Div(bOvO + rvbGvG);
else
    oil = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) + s.Div(bOvO);
end

% Conservation of mass for gas. Again, we have two cases depending on
% whether the model allows us to dissolve the gas phase into the oil phase.
if model.disgas
    % The gas transported in the oil phase.
    rsbOvO = s.faceUpstr(upco, rs).*bOvO;

    gas = (s.pv/dt).*( pvMult.* (bG.* sG  + rs.* bO.* sO) - ...
        pvMult0.*(bG0.*sG0 + rs0.*bO0.*sO0 ) ) + ...
        s.Div(bGvG + rsbOvO);
else
    gas = (s.pv/dt).*( pvMult.*bG.*sG - pvMult0.*bG0.*sG0 ) + s.Div(bGvG);
end

% polymer in water equation :
poro =  s.pv./G.cells.volumes;
polymer = (s.pv.*(1-f.dps)/dt).*(pvMult.*bW.*sW.*cp - ...
   pvMult0.*f.bW(p0).*sW0.*cp0) + (s.pv/dt).* ...
   ( f.rhoR.*((1-poro)./poro).*(ads - ads0)) + s.Div(bWvP);

% Applying correction to the polymer equation when the Jacobian is
% prolematic for some cells.
% Typically it is due to totally and almost non-existence of water.
if ~opt.resOnly
    epsilon = 1.e-8;
    % the first way is based on the diagonal values of the resulting
    % Jacobian matrix
    eps = sqrt(epsilon)*mean(abs(diag(polymer.jac{4})));
    % sometimes there is no water in the whole domain
    if (eps == 0.)
        eps = epsilon;
    end
    % bad marks the cells prolematic in evaluating Jacobian
    bad = abs(diag(polymer.jac{4})) < eps;
    % the other way is to choose based on the water saturation
    polymer(bad) = cp(bad);
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
                                                 dissolved, {cp}, ...
                                                 drivingForces);

% Finally, add in and setup well equations
wc    = vertcat(W.cells);
perf2well = getPerforationToWellMapping(W);
sgn = vertcat(W.sign);

wc_inj = wc(sgn(perf2well) > 0);
cpw     = cp(wc_inj);

% remove the old viscosity and applying the fully mixed viscosity
muWMultW = muWMult(wc_inj);
muWFullyMixed = model.fluid.muWMult(cpw);

mob{1}(wc_inj) = mob{1}(wc_inj) ./ muWFullyMixed .* muWMultW;


if model.usingShear || model.usingShearLog || model.usingShearLogshrate
    % applying the shear effects
    mob{1}(wc) = mob{1}(wc)./shearMultW;
end

[eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, types, wellSol0, wellSol, wellVars, wellMap, p, mob, rho, dissolved, {cp}, dt, opt);
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

end

%--------------------------------------------------------------------------

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

%--------------------------------------------------------------------------

% Effective adsorption, depending of desorption or not
function y = effads(cp, cpmax, f)
   if f.adsInx == 2
      y = f.ads(max(cp, cpmax));
   else
      y = f.ads(cp);
   end
end

%--------------------------------------------------------------------------
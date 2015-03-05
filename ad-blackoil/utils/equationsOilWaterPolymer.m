function [problem, state] = equationsOilWaterPolymer(state0, state, ...
   model, dt, drivingForces, varargin)
% Get linearized problem for oil/water/polymer system with black oil-style
% properties
opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'iteration', -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.Wells;
s = model.operators;

% Properties at current timestep
[p, sW, c, cmax, wellSol] = model.getProps(state, 'pressure', 'water', ...
    'polymer', 'polymermax', 'wellsol');

% Properties at previous timestep
[p0, sW0, c0, cmax0] = model.getProps(state0, 'pressure', 'water', ...
   'polymer', 'polymermax');

pBH    = vertcat(wellSol.bhp);
qWs    = vertcat(wellSol.qWs);
qOs    = vertcat(wellSol.qOs);
qWPoly = vertcat(wellSol.qWPoly);

% Initialize independent variables.
if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [p, sW, c, qWs, qOs, qWPoly, pBH] = ...
            initVariablesADI(p, sW, c, qWs, qOs, qWPoly, pBH);
    else
        zw = zeros(size(pBH));
        [p0, sW0, c0, zw, zw, zw, zw] = ...
            initVariablesADI(p0, sW0, c0, zw, zw, zw, zw); %#ok
        clear zw
    end
end

% We will solve for pressure, water saturation (oil saturation follows via
% the definition of saturations), polymer concentration and well rates +
% bhp.
primaryVars = {'pressure', 'sW', 'polymer', 'qWs', 'qOs', 'qWPoly', 'bhp'};

% Evaluate relative permeability
sO  = 1 - sW;
sO0 = 1 - sW0;

[krW, krO] = model.evaluteRelPerm({sW, sO});

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

% Modifiy relperm by mobility multiplier (if any)
krW = mobMult.*krW; krO = mobMult.*krO;

% Compute transmissibility
T = s.T.*transMult;

% Gravity contribution
gdz = model.getGravityGradient();

% Evaluate water and polymer properties
ads  = effads(c, cmax, model);
ads0 = effads(c0, cmax0, model);
[vW, vP, bW, ~, mobW, mobP, rhoW, pW, upcw, a] = ...
    getFluxAndPropsWaterPolymer_BO(model, p, sW, c, ads, ...
    krW, T, gdz);
bW0 = model.fluid.bW(p0);

% Evaluate oil properties
[vO, bO, mobO, rhoO, p, upco] = getFluxAndPropsOil_BO(model, p, sO, krO, T, gdz);
bO0 = getbO_BO(model, p0);

if model.outputFluxes
    state = model.storeFluxes(state, vW, vO, vP);
end
if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, []);
    state = model.storeMobilities(state, mobW, mobO, mobP);
    state = model.storeUpstreamIndices(state, upcw, upco, []);
end

% EQUATIONS ---------------------------------------------------------------
% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.
bOvO = s.faceUpstr(upco, bO).*vO;
bWvW = s.faceUpstr(upcw, bW).*vW;
bWvP = s.faceUpstr(upcw, bW).*vP;

% Conservation of mass for water
water = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + s.Div(bWvW);

% Conservation of mass for oil
oil = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) + s.Div(bOvO);

% Conservation of polymer in water:
poro = model.rock.poro;
f    = model.fluid;
polymer = (s.pv.*(1-f.dps)/dt).*(pvMult.*bW.*sW.*c - ...
   pvMult0.*bW0.*sW0.*c0) + (s.pv/dt).* ...
   ( f.rhoR.*((1-poro)./poro).*(ads-ads0) ) + s.Div(bWvP);

eqs   = {water, oil, polymer};
names = {'water', 'oil', 'polymer'};
types = {'cell', 'cell', 'cell'};

% Add in any fluxes / source terms prescribed as boundary conditions.
[eqs, qBC, BCTocellMap, qSRC, srcCells] = addFluxesFromSourcesAndBC(...
   model, eqs, {pW, p}, {rhoW, rhoO}, {mobW, mobO}, {bW, bO},  ...
   {sW, sO}, drivingForces);

% Add polymer boundary conditions
if ~isempty(drivingForces.bc) && isfield(drivingForces.bc, 'poly')
   injInx  = qBC{1} > 0; % water inflow indecies
   cbc     = (BCTocellMap')*c; % BCTocellMap' = cellToBCMap
   cbc(injInx) = drivingForces.bc.poly(injInx);
   eqs{3}  = eqs{3} - BCTocellMap*(cbc.*qBC{1});
end

% Add polymer source
if ~isempty(drivingForces.src) && isfield(drivingForces.src, 'poly')
   injInx  = qSRC{1} > 0; % water inflow indecies
   csrc    = c(srcCells);
   csrc(injInx) = drivingForces.src.poly(injInx);
   eqs{3}(srcCells) = eqs{3}(srcCells) - csrc.*qSRC{1};
end

% Finally, add in and setup well equations
if ~isempty(W)
    wm = model.wellmodel;
    if ~opt.reverseMode
        wc   = vertcat(W.cells);
        pw   = p(wc);
        rhos = [f.rhoWS, f.rhoOS];
        bw   = {bW(wc), bO(wc)};
        mw   = {mobW(wc), mobO(wc)};
        s    = {sW(wc), sO(wc)};

        [cqs, weqs, ctrleqs, wc, state.wellSol] = ...
            wm.computeWellFlux(model, W, wellSol, ...
            pBH, {qWs, qOs}, pw, rhos, bw, mw, s, {},...
            'nonlinearIteration', opt.iteration);

        % Store the well equations (relate well bottom hole pressures to
        % influx).
        eqs(4:5) = weqs;
        % Store the control equations (trivial equations ensuring that each
        % well will have values corresponding to the prescribed value)
        eqs{7} = ctrleqs;
        % Add source terms to the equations. Negative sign may be
        % surprising if one is used to source terms on the right hand side,
        % but this is the equations on residual form.
        eqs{1}(wc) = eqs{1}(wc) - cqs{1};
        eqs{2}(wc) = eqs{2}(wc) - cqs{2};

        % Polymer well equations
        [~, wciPoly, iInxW] = getWellPolymer(W);
        cw        = c(wc);
        cw(iInxW) = wciPoly;
        cbarw     = cw/f.cmax;

        % Divide away water mobility and add in polymer
        bWqP = cw.*cqs{1}./(a + (1-a).*cbarw);
        eqs{3}(wc) = eqs{3}(wc) - bWqP;

        % Well polymer rate for each well is water rate in each perforation
        % multiplied with polymer concentration in that perforated cell.
        perf2well = getPerforationToWellMapping(W);
        Rw = sparse(perf2well, (1:numel(perf2well))', 1, ...
           numel(W), numel(perf2well));
        eqs{6} = qWPoly - Rw*(cqs{1}.*cw);

        names(4:7) = {'waterWells', 'oilWells', 'polymerWells', ...
            'closureWells'};
        types(4:7) = {'perf', 'perf', 'perf', 'well'};
    else
        [eq, n, typ] = ...
            wm.createReverseModeWellEquations(model, state0.wellSol, p0);
        % Add another equation for polymer well rates
        [eqs{4:7}] = deal(eq{1});
        [names{4:7}] = deal(n{1});
        [types{4:7}] = deal(typ{1});
    end
end
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
    polInj = cellfun(@(x)~isempty(x), {W(inj).poly});
    wPoly = zeros(nnz(inj), 1);
    wPoly(polInj) = vertcat(W(inj(polInj)).poly);
    wciPoly = rldecode(wPoly, cellfun(@numel, {W(inj).cells}));

    % Injection cells
    nPerf = cellfun(@numel, {W.cells})';
    nw    = numel(W);
    perf2well = rldecode((1:nw)', nPerf);
    compi = vertcat(W.compi);
    iInx  = rldecode(inj, nPerf);
    iInx  = find(iInx);
    iInxW = iInx(compi(perf2well(iInx),1)==1);
end

% Effective adsorption, depending of desorption or not
function y = effads(c, cmax, model)
   if model.fluid.adsInx == 2
      y = model.fluid.ads(max(c, cmax));
   else
      y = model.fluid.ads(c);
   end
end




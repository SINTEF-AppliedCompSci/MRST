function [problem, state] = pressureEquationOilWaterPolymer(state0, ...
    state, model, dt, drivingForces, varargin)

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'propsPressure', [], ...
             'iteration', -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.Wells;
perf2well = getPerforationToWellMapping(W);

assert(isempty(drivingForces.bc) && isempty(drivingForces.src))

s = model.operators;
f = model.fluid;

[p, sW, c, cmax, wellSol] = model.getProps(state, 'pressure', 'water', ...
    'polymer', 'polymermax', 'wellsol');
[p0, sW0, c0, cmax0] = model.getProps(state0, 'pressure', 'water', ...
   'polymer', 'polymermax');


pBH    = vertcat(wellSol.bhp);
qWs    = vertcat(wellSol.qWs);
qOs    = vertcat(wellSol.qOs);
qWPoly = vertcat(wellSol.qWPoly);

%Initialization of independent variables ----------------------------------

if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [p, qWs, qOs, qWPoly, pBH] = ...
            initVariablesADI(p, qWs, qOs, qWPoly, pBH);
    else
        assert(0, 'Backwards solver not supported for splitting');
    end
end
primaryVars = {'pressure', 'qWs', 'qOs', 'qWPoly', 'bhp'};

p_prop = opt.propsPressure;
otherPropPressure = ~isempty(p_prop);
if ~otherPropPressure
    p_prop = p;
end

% -------------------------------------------------------------------------
sO  = 1 - sW;
sO0 = 1 - sW0;

[krW, krO] = model.evaluteRelPerm({sW, sO});

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, ...
    p_prop, p0);

% Modifiy relperm by mobility multiplier (if any)
krW = mobMult.*krW; krO = mobMult.*krO;

% Compute transmissibility
T = s.T.*transMult;

% Gravity contribution
gdz = model.getGravityGradient();

% Evaluate water and polymer properties
ads  = effads(c, cmax, model);
[vW, vP, bW, muWMult, mobW, mobP, rhoW, pW, upcw, a, dpW] = ...
    getFluxAndPropsWaterPolymer_BO(model, p_prop, sW, c, ads, ...
    krW, T, gdz);
bW0 = model.fluid.bW(p0);

% Evaluate oil properties
[vO, bO, mobO, rhoO, pO, upco, dpO] = getFluxAndPropsOil_BO(model, ...
    p_prop, sO, krO, T, gdz);
bO0 = getbO_BO(model, p0);

if otherPropPressure
    % We have used a different pressure for property evaluation, undo the
    % effects of this on the fluxes.
    dp_diff = s.Grad(p) - s.Grad(p_prop);
    
    vW = -s.faceUpstr(upcw, mobW).*s.T.*(dpW + dp_diff);
    vO = -s.faceUpstr(upco, mobO).*s.T.*(dpO + dp_diff);
    vP = -s.faceUpstr(upcw, mobP).*s.T.*(dpW + dp_diff);
end

% These are needed in transport solver, so we output them regardless of
% any flags set in the model.
state = model.storeFluxes(state, vW, vO, vP);
state = model.storeUpstreamIndices(state, upcw, upco, upcw);
if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, bW);
    state = model.storeMobilities(state, mobW, mobO, mobP);
end

% EQUATIONS ---------------------------------------------------------------
% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.

bOvO = s.faceUpstr(upco, bO).*vO;
bWvW = s.faceUpstr(upcw, bW).*vW;

% Oil
oil = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0) + s.Div(bOvO);

% Water
wat = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + s.Div(bWvW);

[eqs, names, types] = deal({});

% well equations
if ~isempty(W)
    wc   = vertcat(W.cells);
    pw   = p(wc);
    rhos = [f.rhoWS, f.rhoOS];
    bw   = {bW(wc), bO(wc)};
    mw   = {mobW(wc), mobO(wc)};
    s    = {sW(wc), 1-sW(wc)};

    wm = model.wellmodel;
    [cqs, weqs, ctrleqs, wc, state.wellSol, cqr] = ...
        wm.computeWellFlux(model, W, wellSol, pBH, {qWs, qOs}, pw, ...
        rhos, bw, mw, s, {}, 'nonlinearIteration', opt.iteration);
    eqs(2:3) = weqs;
    eqs{5} = ctrleqs;
    
    qW = cqr{1};
    qO = cqr{2};
    
    wat(wc) = wat(wc) - cqs{1};
    oil(wc) = oil(wc) - cqs{2};
    
    % Polymer well equations
    % Well polymer rate for each well is water rate in each perforation
    % multiplied with polymer concentration in that perforated cell.
    [~, wciPoly, iInxW] = getWellPolymer(W);
    cw        = c(wc);
    cw(iInxW) = wciPoly;
    Rw = sparse(perf2well, (1:numel(perf2well))', 1, ...
       numel(W), numel(perf2well));
    eqs{4} = qWPoly - Rw*(cqs{1}.*cw);
    
    names(2:5) = {'oilWells', 'waterWells', 'polymerWells', ...
        'closureWells'};
    types(2:5) = {'perf', 'perf', 'perf', 'well'};
end

eqs{1}   = oil./bO + wat./bW;
names{1} = 'pressure';
types{1} = 'cell';

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

for i = 1:numel(W)
    wp = perf2well == i;
    state.wellSol(i).flux = [double(qW(wp)), double(qO(wp))];
end

state.s0 = state0.s;
state.c0 = state0.c;
state.bfactor0 = [double(bW0), double(bO0)];

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




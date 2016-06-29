function [problem, state] = pressureEquationOilWater(state0, state, model, dt, drivingForces, varargin)

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'propsPressure', [], ...
             'staticWells',  false, ...
             'iteration', -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.W;

% assert(isempty(drivingForces.bc) && isempty(drivingForces.src))

s = model.operators;
f = model.fluid;
G = model.G;

[p, sW, wellSol] = model.getProps(state, 'pressure', 'water', 'wellsol');
[p0, sW0] = model.getProps(state0, 'pressure', 'water');


pBH    = vertcat(wellSol.bhp);
qWs    = vertcat(wellSol.qWs);
qOs    = vertcat(wellSol.qOs);

%Initialization of independent variables ----------------------------------

if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [p, qWs, qOs, pBH] = ...
            initVariablesADI(p, qWs, qOs, pBH);
    else
        assert(0, 'Backwards solver not supported for splitting');
    end
end
primaryVars = {'pressure', 'qWs', 'qOs', 'bhp'};

p_prop = opt.propsPressure;
otherPropPressure = ~isempty(p_prop);
if ~otherPropPressure
    p_prop = p;
end

% -------------------------------------------------------------------------
sO  = 1 - sW;
sO0 = 1 - sW0;

[krW, krO] = model.evaluateRelPerm({sW, sO});

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p_prop, p0);

% Modifiy relperm by mobility multiplier (if any)
krW = mobMult.*krW; krO = mobMult.*krO;

% Compute transmissibility
T = s.T.*transMult;

% Gravity contribution
gdz = model.getGravityGradient();


% Evaluate water properties
[vW, bW, mobW, rhoW, pW, upcw, dpW] = getFluxAndPropsWater_BO(model, p_prop, sW, krW, T, gdz);
bW0 = f.bW(p0);

% Evaluate oil properties
[vO, bO, mobO, rhoO, pO, upco, dpO] = getFluxAndPropsOil_BO(model, p_prop, sO, krO, T, gdz);
bO0 = getbO_BO(model, p0);

if otherPropPressure
    % We have used a different pressure for property evaluation, undo the
    % effects of this on the fluxes.
    dp_diff = s.Grad(p) - s.Grad(p_prop);
    
    vW = -s.faceUpstr(upcw, mobW).*s.T.*(dpW + dp_diff);
    vO = -s.faceUpstr(upco, mobO).*s.T.*(dpO + dp_diff);
end

% These are needed in transport solver, so we output them regardless of
% any flags set in the model.
state = model.storeFluxes(state, vW, vO, []);
state = model.storeUpstreamIndices(state, upcw, upco, []);
if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, []);
    state = model.storeMobilities(state, mobW, mobO, []);
end
% EQUATIONS ---------------------------------------------------------------
% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.
bOvO = s.faceUpstr(upco, bO).*vO;
bWvW = s.faceUpstr(upcw, bW).*vW;


oil = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0) + s.Div(bOvO);
% water:
wat = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + s.Div(bWvW);

eqTmp = {wat, oil};
[eqTmp, ~, qRes] = addFluxesFromSourcesAndBC(model, eqTmp, ...
                                       {pW, p},...
                                       {rhoW,     rhoO},...
                                       {mobW,     mobO}, ...
                                       {bW, bO},  ...
                                       {sW, sO}, ...
                                       drivingForces);
wat = eqTmp{1};
oil = eqTmp{2};

if model.outputFluxes
    state = model.storeBoundaryFluxes(state, qRes{1}, qRes{2}, [], drivingForces);
end
[eqs, names, types] = deal({});

% well equations
if ~isempty(W)
    wc    = vertcat(W.cells);
    perf2well = getPerforationToWellMapping(W);
    if opt.staticWells
        q = vertcat(state.wellSol.flux);
        
        qW = q(:, 1);
        qO = q(:, 2);
        
        cqs = {bW(wc).*qW, bO(wc).*qO};
    else
        pw   = p(wc);
        rhos = [f.rhoWS, f.rhoOS];
        bw   = {bW(wc), bO(wc)};
        mw   = {mobW(wc), mobO(wc)};
        sat = {sW(wc), 1 - sW(wc)};

        wm = model.wellmodel;
        [cqs, weqs, ctrleqs, wc, state.wellSol, cqr]  = wm.computeWellFlux(model, W, wellSol, ...
                                             pBH, {qWs, qOs}, pw, rhos, bw, mw, sat, {},...
                                             'nonlinearIteration', opt.iteration);
        eqs(2:3) = weqs;
        eqs{4} = ctrleqs;

        qW = cqr{1};
        qO = cqr{2};
        
        names(2:4) = {'oilWells', 'waterWells', 'closureWells'};
        types(2:4) = {'perf', 'perf', 'well'};

    end
    
    wat(wc) = wat(wc) - cqs{1};
    oil(wc) = oil(wc) - cqs{2};
end

eqs{1} = (dt./s.pv).*(oil./bO + wat./bW);
names{1} = 'pressure';
types{1} = 'cell';

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

for i = 1:numel(W)
    wp = perf2well == i;
    state.wellSol(i).flux = [double(qW(wp)), double(qO(wp))];
end

end

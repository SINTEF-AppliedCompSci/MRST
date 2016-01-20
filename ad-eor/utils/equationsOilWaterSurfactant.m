function [problem, state] = equationsOilWaterSurfactant(state0, state, ...
   model, dt, drivingForces, varargin)
% Get linearized problem for oil/water/surfactant system with black oil-style
% properties
opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'iteration', -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.W;
s = model.operators;
fluid = model.fluid;

% Properties at current timestep
[p, sW, c, wellSol] = model.getProps(state, 'pressure', 'water', 'surfactant', 'wellsol');

% Properties at previous timestep
[p0, sW0, c0] = model.getProps(state0, 'pressure', 'water', 'surfactant');

pBH    = vertcat(wellSol.bhp);
qWs    = vertcat(wellSol.qWs);
qOs    = vertcat(wellSol.qOs);
qWSft = vertcat(wellSol.qWSft);

% Initialize independent variables.
if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [p, sW, c, qWs, qOs, qWSft, pBH] = ...
            initVariablesADI(p, sW, c, qWs, qOs, qWSft, pBH);
    else
        zw = zeros(size(pBH));
        [p0, sW0, c0, zw, zw, zw, zw] = ...
            initVariablesADI(p0, sW0, c0, zw, zw, zw, zw); %#ok
        clear zw
    end
end

% We will solve for pressure, water saturation (oil saturation follows via
% the definition of saturations), surfactant concentration and well rates +
% bhp.
primaryVars = {'pressure', 'sW', 'surfactant', 'qWs', 'qOs', 'qWSft', 'bhp'};

% Evaluate relative permeability
sO  = 1 - sW;
sO0 = 1 - sW0;

Nc = computeCapillaryNumber(p, c, fluid, s);
[krW, krO] = computeRelPermSft(sW, Nc, fluid);

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

% Modifiy relperm by mobility multiplier (if any)
krW = mobMult.*krW; 
krO = mobMult.*krO;

% Compute transmissibility
T = s.T.*transMult;

% Gravity contribution
gdz = model.getGravityGradient();

% NOT IMPLEMENTED YET: 
% Adsortion.
% Viscosity change.
% Capillary pressue
pW = p;
pO = p;

fluid = model.fluid;

bW0 = fluid.bW(p0);
bO0 = fluid.bO(p0);


% Water and Surfactant flux
bW   = fluid.bW(p);
rhoW = bW.*fluid.rhoWS;
rhoWf  = s.faceAvg(rhoW);
muW  = fluid.muW(p);
mobW = krW./muW;
dpW  = s.Grad(pW) - rhoWf.*gdz;
upcw = (double(dpW)<=0);
vW   = -s.faceUpstr(upcw, mobW).*s.T.*dpW;
mobSft = mobW.*c;
vSft   = - s.faceUpstr(upcw, mobSft).*s.T.*dpW;

% Oil flux
bO   = fluid.bO(pO);
rhoO = bO.*fluid.rhoOS;
rhoOf  = s.faceAvg(rhoO);
if isfield(fluid, 'BOxmuO')
   muO = fluid.BOxmuO(pO).*bO;
else
   muO = fluid.muO(pO);
end
mobO = krO./muO;
dpO  = s.Grad(pO) - rhoOf.*gdz;
upco = (double(dpO)<=0);
vO   = -s.faceUpstr(upco, mobO).*s.T.*dpO;


if model.outputFluxes 
    state = model.storeFluxes(state, vW, vO, vSft);
end
if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, []);
    state = model.storeMobilities(state, mobW, mobO, mobSft);
    state = model.storeUpstreamIndices(state, upcw, upco, []);
end

% EQUATIONS ---------------------------------------------------------------
% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.
bOvO = s.faceUpstr(upco, bO).*vO;
bWvW = s.faceUpstr(upcw, bW).*vW;
bWvSft = s.faceUpstr(upcw, bW).*vSft;

% Conservation of mass for water
water = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + s.Div(bWvW);

% Conservation of mass for oil
oil = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) + s.Div(bOvO);

% Conservation of surfactant in water:
poro = model.rock.poro;
surfactant = (s.pv/dt).*(pvMult.*bW.*sW.*c - pvMult0.*bW0.*sW0.*c0) + s.Div(bWvSft);

eqs   = {water, oil, surfactant};
names = {'water', 'oil', 'surfactant'};
types = {'cell', 'cell', 'cell'};

% Add in any fluxes / source terms prescribed as boundary conditions.
[eqs, qBC, BCTocellMap, qSRC, srcCells] = addFluxesFromSourcesAndBC(...
   model, eqs, {pW, p}, {rhoW, rhoO}, {mobW, mobO}, {bW, bO},  ...
   {sW, sO}, drivingForces);

% Add surfactant boundary conditions
% setup functions not implemented
if ~isempty(drivingForces.bc) && isfield(drivingForces.bc, 'surfact')
   injInx  = qBC{1} > 0; % water inflow indecies
   cbc     = (BCTocellMap')*c; % BCTocellMap' = cellToBCMap
   cbc(injInx) = drivingForces.bc.surfact(injInx); 
   eqs{3}  = eqs{3} - BCTocellMap*(cbc.*qBC{1});
end

% Add surfactant source
% setup functions not implemented 
if ~isempty(drivingForces.src) && isfield(drivingForces.src, 'surfact')
   injInx  = qSRC{1} > 0; % water inflow indecies
   csrc    = c(srcCells);
   csrc(injInx) = drivingForces.src.surfact(injInx);
   eqs{3}(srcCells) = eqs{3}(srcCells) - csrc.*qSRC{1};
end

% Finally, add in and setup well equations
if ~isempty(W)
    wm = model.wellmodel;
    if ~opt.reverseMode
        wc   = vertcat(W.cells);
        pw   = p(wc);
        rhos = [fluid.rhoWS, fluid.rhoOS];
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

        % surfactant well equations
        [~, wciSft, iInxW] = getWellSurfactant(W);
        cw        = c(wc);
        cw(iInxW) = wciSft;

        % Add surfactant
        bWqP = cw.*cqs{1};
        eqs{3}(wc) = eqs{3}(wc) - bWqP;

        % Well surfactant rate for each well is water rate in each perforation
        % multiplied with surfactant concentration in that perforated cell.
        perf2well = getPerforationToWellMapping(W);
        Rw = sparse(perf2well, (1:numel(perf2well))', 1, ...
           numel(W), numel(perf2well));
        eqs{6} = qWSft - Rw*(cqs{1}.*cw);

        names(4:7) = {'waterWells', 'oilWells', 'surfactantWells', 'closureWells'};
        types(4:7) = {'perf', 'perf', 'perf', 'well'};
    else
        [eq, n, typ] = ...
            wm.createReverseModeWellEquations(model, state0.wellSol, p0);
        % Add another equation for surfactant well rates
        [eqs{4:7}] = deal(eq{1});
        [names{4:7}] = deal(n{1});
        [types{4:7}] = deal(typ{1});
    end
end
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end


%--------------------------------------------------------------------------


function [wSft, wciSft, iInxW] = getWellSurfactant(W)
    if isempty(W)
        wSft = [];
        wciSft = [];
        iInxW = [];
        return
    end
    inj   = vertcat(W.sign)==1;
    surfactInj = cellfun(@(x)~isempty(x), {W(inj).surfact});
    wSft = zeros(nnz(inj), 1);
    wSft(surfactInj) = vertcat(W(inj(surfactInj)).surfact);
    wciSft = rldecode(wSft, cellfun(@numel, {W(inj).cells}));

    % Injection cells
    nPerf = cellfun(@numel, {W.cells})';
    nw    = numel(W);
    perf2well = rldecode((1:nw)', nPerf);
    compi = vertcat(W.compi);
    iInx  = rldecode(inj, nPerf);
    iInx  = find(iInx);
    iInxW = iInx(compi(perf2well(iInx),1)==1);
end





function [problem, state] = equationsWater(state0, state, model, dt, drivingForces, varargin)
% Get linearized problem for single phase water system with black oil-style
% properties

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'iteration', -1);  % Compatibility only

opt = merge_options(opt, varargin{:});

W = drivingForces.W;
%assert(isempty(drivingForces.bc) && isempty(drivingForces.src))
assert(isempty(drivingForces.src))


s = model.operators;
G = model.G;
f = model.fluid;

[p, wellSol] = model.getProps(state, 'pressure', 'wellsol');

[p0] = model.getProps(state0, 'pressure');

pBH    = vertcat(wellSol.bhp);
qWs    = vertcat(wellSol.qWs);

%Initialization of independent variables ----------------------------------

if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [p, qWs, pBH] = ...
            initVariablesADI(p, qWs, pBH);
    else
        [p0, tmp, tmp] = ...
            initVariablesADI(p0, sW0,          ...
            zeros(size(qWs)), ...
            zeros(size(qOs)), ...
            zeros(size(pBH)));                          %#ok
    end
end
primaryVars = {'pressure','qWs', 'bhp'};

clear tmp
%grav  = gravity;
%gdz   = s.Grad(G.cells.centroids) * model.gravity';
gdz   = s.Grad(G.cells.centroids) * model.getGravityVector()';
%--------------------
%check for p-dependent tran mult:
trMult = 1;
if isfield(f, 'tranMultR'), trMult = f.tranMultR(p); end

%check for p-dependent porv mult:
pvMult = 1; pvMult0 = 1;
if isfield(f, 'pvMultR')
    pvMult =  f.pvMultR(p);
    pvMult0 = f.pvMultR(p0);
end
transMult=1;
if isfield(f, 'transMult')
   transMult=f.transMult(p); 
end

trans=s.T.*transMult;
% -------------------------------------------------------------------------
% water props (calculated at oil pressure OK?)
bW     = f.bW(p);
rhoW   = bW.*f.rhoWS;
% rhoW on face, avarge of neighboring cells (E100, not E300)
rhoWf  = s.faceAvg(rhoW);
mobW   = trMult./f.muW(p);
dpW     = s.Grad(p) - rhoWf.*gdz;
% water upstream-index
upcw = (double(dpW)<=0);
vW = - s.faceUpstr(upcw, mobW).*trans.*dpW;
bWvW = s.faceUpstr(upcw, bW).*vW;

if model.outputFluxes
    state = model.storeFluxes(state, vW, [], []);
end

if model.extraStateOutput
    state = model.storebfactors(state, bW, [], []);
    state = model.storeMobilities(state, mobW, [], []);
    state = model.storeUpstreamIndices(state, upcw, [], []);
end
% EQUATIONS ---------------------------------------------------------------
% water:
eqs{1} = (s.pv/dt).*( pvMult.*bW - pvMult0.*f.bW(p0) ) + s.Div(bWvW);
eqs = addFluxesFromSourcesAndBC(model, eqs, ...
                                       {p},...
                                       {rhoW},...
                                       {mobW}, ...
                                       {bW},  ...
                                       {1}, ...
                                       drivingForces);

names = {'water'};
types = {'cell'};
% well equations
if ~isempty(W)
    if ~opt.reverseMode
        wc    = vertcat(W.cells);
        pw   = p(wc);
        rhos = [f.rhoWS];
        bw   = {bW(wc)};
        mw   = {mobW(wc)};
        s = {1};

        wm = WellModel();
        [cqs, weqs, ctrleqs, wc, state.wellSol]  = wm.computeWellFlux(model, W, wellSol, ...
                                             pBH, {qWs}, pw, rhos, bw, mw, s, {},...
                                             'nonlinearIteration', opt.iteration,'referencePressureIndex', 1);
        eqs(2) = weqs;
        eqs{3} = ctrleqs;
        
        eqs{1}(wc) = eqs{1}(wc) - cqs{1};        
        names(2:3) = {'waterWells', 'closureWells'};
        types(2:3) = {'perf', 'well'};
    else
        % in reverse mode just gather zero-eqs of correct size
        for eqn = 2:3
            nw = numel(state0.wellSol);
            zw = double2ADI(zeros(nw,1), p0);
            eqs(2:3) = {zw, zw};
        end
        names(2:3) = {'empty', 'empty'};
        types(2:3) = {'none', 'none'};
    end
end
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end
%--------------------------------------------------------------------------










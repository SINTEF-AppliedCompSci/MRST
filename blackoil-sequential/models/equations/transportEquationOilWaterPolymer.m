function [problem, state] = transportEquationOilWaterPolymer(state0, state, model, dt, drivingForces, varargin)

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'scaling', [],...
             'resOnly', false,...
             'history', [],...
             'solveForWater', false, ...
             'solveForOil', true, ...
             'iteration', -1, ...
             'stepOptions', []);  % Compatibility only

opt = merge_options(opt, varargin{:});


W = drivingForces.Wells;
assert(isempty(drivingForces.bc) && isempty(drivingForces.src))

s = model.operators;
f = model.fluid;
G = model.G;

assert(~(opt.solveForWater && opt.solveForOil));

% Properties at current timestep
[p, sW, c, cmax, wellSol] = model.getProps(state, 'pressure', 'water', ...
    'polymer', 'polymermax', 'wellsol');

% Properties at previous timestep
[p0, sW0, c0, cmax0] = model.getProps(state0, 'pressure', 'water', ...
   'polymer', 'polymermax');

wflux = sum(vertcat(wellSol.flux), 2);

%Initialization of independent variables ----------------------------------

if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [sW, c] = initVariablesADI(sW, c);
    else
        assert(0, 'Backwards solver not supported for splitting');
    end
end
primaryVars = {'sW', 'polymer'};

clear tmp

% -------------------------------------------------------------------------

sO = 1 - sW;

[krW, krO] = model.evaluteRelPerm({sW, sO});

% Gravity contribution
gdz = model.getGravityGradient();

% Evaluate water and polymer properties
ads  = effads(c, cmax, model);
ads0 = effads(c0, cmax0, model);
[~, ~, bW, ~, mobW, mobP, rhoW, pW, upcw, a] = ...
    getFluxAndPropsWaterPolymer_BO(model, p, sW, c, ads, ...
    krW, s.T, gdz);

% Water
%[bW, rhoW, mobW, Gw] = propsOW_water(sW, krW, gdz, f, p, s);

% Oil properties
[bO, rhoO, mobO, Go] = propsOW_oil(  sO, krO, gdz, f, p, s);

    
if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, []);
    state = model.storeMobilities(state, mobW, mobO, mobP);
end

if ~isempty(W)
    perf2well = getPerforationToWellMapping(W);
    wc = vertcat(W.cells);
    
    mobWw = mobW(wc);
    mobOw = mobO(wc);
    totMobw = mobWw + mobOw;

    f_w_w = mobWw./totMobw;
    f_o_w = mobOw./totMobw;

    isInj = wflux > 0;
    compWell = vertcat(W.compi);
    compPerf = compWell(perf2well, :);

    f_w_w(isInj) = compPerf(isInj, 1);
    f_o_w(isInj) = compPerf(isInj, 2);
    
    bWqW = bW(wc).*f_w_w.*wflux;
    bOqO = bO(wc).*f_o_w.*wflux;
    
    % Polymer well equations
    [~, wciPoly, iInxW] = getWellPolymer(W);
    cw        = c(wc);
    cw(iInxW) = wciPoly;
    cbarw     = cw/f.cmax;
    
    % Divide away water mobility and add in polymer
    mcw  = cw./(a + (1-a).*cbarw);
	bWqP = bWqW.*mcw;

    % Store well fluxes
    wflux_O = double(bOqO);
    wflux_W = double(bWqW);
    wflux_P = double(bWqP);
    
    for i = 1:numel(W)
        perfind = perf2well == i;
        state.wellSol(i).qOs = sum(wflux_O(perfind));
        state.wellSol(i).qWs = sum(wflux_W(perfind));
        state.wellSol(i).qPs = sum(wflux_P(perfind)); % ??
    end

end
%check for p-dependent porv mult:
pvMult = 1; pvMult0 = 1;
if isfield(f, 'pvMultR')
    pvMult =  f.pvMultR(p);
    pvMult0 = f.pvMultR(p0);
end


% Get total flux from state
flux = sum(state.flux, 2);
vT = flux(model.operators.internalConn);

% Stored upstream indices
% if model.staticUpwind
%     flag = state.upstreamFlag;
% else
%     flag = multiphaseUpwindIndices({Gw, Go}, vT, s.T, {mobW, mobO}, s.faceUpstr);
% end
% 
% upcw  = flag(:, 1);
% upco  = flag(:, 2);

% oil upstream-index
upco = (double(dpO)<=0);

mobOf = s.faceUpstr(upco, mobO);
mobWf = s.faceUpstr(upcw, mobW);

% m(c) = muweff / mupeff => vwp = m(c)*vw
mc = c./(a + (1-a).*(c/f.cmax));

totMob = (mobOf + mobWf);
totMob = max(totMob, sqrt(eps));

if opt.solveForWater
    f_w = mobWf./totMob;
    bWvW = s.faceUpstr(upcw, bW).*f_w.*(vT + s.T.*mobOf.*(Gw - Go));
    bWvP = bWvW.*mc;
    
    wat = (s.pv/dt).*(pvMult.*bW.*sW       - pvMult0.*f.bW(p0).*sW0    ) + s.Div(bWvW);
    wat(wc) = wat(wc) - bWqW;
    
    % Conservation of polymer in water:
    poro = model.rock.poro;
    f    = model.fluid;
    poly = (s.pv.*(1-f.dps)/dt).*(pvMult.*bW.*sW.*c - ...
       pvMult0.*f.bW(p0).*sW0.*c0) + (s.pv/dt).* ...
       ( f.rhoR.*((1-poro)./poro).*(ads-ads0) ) + s.Div(bWvP);
    poly(wc) = poly(wc) - bWqP;
    
    eqs{1} = wat;
    eqs{2} = poly;
    names = {'water','polymer'};
    types = {'cell','cell'};
else
    error('Not supported for polymer')
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




function [problem, state] = transportEquationOilWaterPolymer(state0, ...
    state, model, dt, drivingForces, varargin)

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

% Polymer shear thinning/thickening
usingShear = isfield(f, 'plyshearMult');

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

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

% Modifiy relperm by mobility multiplier (if any)
krW = mobMult.*krW; krO = mobMult.*krO;

% Compute transmissibility
T = s.T.*transMult;

% Gravity gradient per face
gdz = model.getGravityGradient();

% Evaluate water and polymer properties
ads  = effads(c, cmax, model);
ads0 = effads(c0, cmax0, model);
[vW, vP, bW, muWMult, mobW, mobP, rhoW, pW, upcw, dpW, extraOutput] = ...
    getFluxAndPropsWaterPolymer_BO(model, p, sW, c, ads, krW, T, gdz);

% Evaluate oil properties
[vO, bO, mobO, rhoO, pO, upco, dpO] = getFluxAndPropsOil_BO(model, p, ...
    sO, krO, T, gdz);

gp = s.Grad(p);
Gw = gp - dpW;
Go = gp - dpO;

if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, []);
    state = model.storeMobilities(state, mobW, mobO, mobP);
end

if model.extraPolymerOutput
%     state = model.storeShearMultiplier(state, shearMult);
%     state = model.storeEffectiveWaterVisc(state, extraOutput.muWeff);
%     state = model.storeEffectivePolymerVisc(state, extraOutput.muPeff);
%     state = model.storePolymerAdsorption(state, ads);
%     state = model.storeRelpermReductionFactor(state, extraOutput.Rk);
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
    
    % Polymer injection
%     cw = c(wc);
%     wpoly = vertcat(W.poly);
%     wpoly = wpoly(perf2well);
%     cw(isInj) = wpoly(isInj);
    
    % Polymer well equations
    [~, wciPoly, iInxW] = getWellPolymer(W);
    cw        = c(wc);
    cw(iInxW) = wciPoly;
    
    bWqP = cw.*bWqW;
    
    if usingShear
        % Compute shear rate multiplier for wells
        % The water velocity is computed using a the reprensentative 
        % radius rR.
        % rR = sqrt(rW * rA)
        % rW is the well bore radius.
        % rA is the equivalent radius of the grid block in which the 
        %    wellis completed.
        
        [~, wciPoly, iInxW] = getWellPolymer(W);
        
        assert(isfield(W, 'rR'), ...
            'The representative radius needs to be suppplied.');

        muWMultW = muWMult(wc);
        % Maybe should also apply this for PRODUCTION wells.
        muWMultW((iInxW(wciPoly==0))) = 1;

        % The following formulations assume that the wells are always
        % in the z direction 
        % IMPROVED HERE LATER
        [~, ~, dz] = cellDims(model.G, wc);
        
        rR = vertcat(W.rR);
        VW0W = double(bWqW)./(model.rock.poro(wc).*rR.*dz*2*pi);
        shearMultW = getPolymerShearMultiplier(model, VW0W, muWMultW);

        % Apply shear velocity multiplier
        bWqW = bWqW.*shearMultW;
        bWqP = bWqP.*shearMultW;
    end
    
    % Store well fluxes
    wflux_O = double(bOqO);
    wflux_W = double(bWqW);
    wflux_P = double(bWqP);
    
    for i = 1:numel(W)
        perfind = perf2well == i;
        state.wellSol(i).qOs = sum(wflux_O(perfind));
        state.wellSol(i).qWs = sum(wflux_W(perfind));
        state.wellSol(i).qPs = sum(wflux_P(perfind));
    end

end

% Get total flux from state
assert(size(state.flux,2)==3, 'Not the flux expected');
flux = sum(state.flux(:,1:2), 2);
vT = flux(model.operators.internalConn);

% Stored upstream indices
if model.staticUpwind
    flag = state.upstreamFlag;
else
    flag = multiphaseUpwindIndices({Gw, Go}, vT, s.T, {mobW, mobO}, ...
        s.faceUpstr);
end

upcw  = flag(:, 1);
upco  = flag(:, 2);

mobOf = s.faceUpstr(upco, mobO);
mobWf = s.faceUpstr(upcw, mobW);
mobPf = s.faceUpstr(upcw, mobP);


%% TEMP TODO TEST
% Change velocitites due to polymer shear thinning / thickening
if usingShear
    poroFace  = s.faceAvg(model.rock.poro);
    faceArea  = model.G.faces.areas(s.internalConn);
    Vw        = vW./(poroFace .* faceArea);
    muWMultf  = s.faceUpstr(upcw, muWMult);
    shearMult = getPolymerShearMultiplier(model, Vw, muWMultf);
    mobWf     = mobWf .* shearMult;
    mobPf     = mobPf .* shearMult;
end
%%

totMob = (mobOf + mobWf);
totMob = max(totMob, sqrt(eps));

if opt.solveForWater
    f_w = mobWf./totMob;
    bWvW   = s.faceUpstr(upcw, bW).*f_w.*(vT + s.T.*mobOf.*(Gw - Go));

    wat = (s.pv/dt).*(pvMult.*bW.*sW - pvMult0.*f.bW(p0).*sW0) + ...
        s.Div(bWvW);
    wat(wc) = wat(wc) - bWqW;
    
    eqs{1} = wat;
    name1  = 'water';
else
    f_o = mobOf./totMob;
    bOvO   = s.faceUpstr(upco, bO).*f_o.*(vT + s.T.*mobWf.*(Go - Gw));
    
    oil = (s.pv/dt).*( pvMult.*bO.*(1-sW) - ...
        pvMult0.*f.bO(p0).*(1-sW0) ) + s.Div(bOvO);
    oil(wc) = oil(wc) - bOqO;
    
    eqs{1} = oil;
    name1  = 'oil';
end

% Polymer equations
f_p = mobPf./totMob;
bWvP   = s.faceUpstr(upcw, bW).*f_p.*(vT + s.T.*mobOf.*(Gw - Go));

poro = model.rock.poro;
poly = (s.pv.*(1-f.dps)/dt).*(pvMult.*bW.*sW.*c - ...
    pvMult0.*f.bW(p0).*sW0.*c0) + (s.pv/dt).* ...
    ( f.rhoR.*((1-poro)./poro).*(ads-ads0) ) + s.Div(bWvP);
poly(wc) = poly(wc) - bWqP;

% Fix for (almost) zero water in the well
if isa(poly, 'ADI')
   epsilon = 1.e-8;
   epsilon = sqrt(epsilon)*mean(abs(diag(poly.jac{2})));
   bad     = abs(diag(poly.jac{2})) < epsilon;
   poly(bad) = c(bad);
end

eqs{2} = poly./f.cmax; % scale with cmax
names  = {name1 'polymer'};
types  = {'cell',  'cell'};

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



function [dx, dy, dz] = cellDims(G, ix)
% cellDims -- Compute physical dimensions of all cells in single well
%
% SYNOPSIS:
%   [dx, dy, dz] = cellDims(G, ix)
%
% PARAMETERS:
%   G  - Grid data structure.
%   ix - Cells for which to compute the physical dimensions
%
% RETURNS:
%   dx, dy, dz -- [dx(k) dy(k)] is bounding box in xy-plane, while dz(k) =
%                 V(k)/dx(k)*dy(k)

    n = numel(ix);
    [dx, dy, dz] = deal(zeros([n, 1]));

    ixc = G.cells.facePos;
    ixf = G.faces.nodePos;

    for k = 1 : n,
       c = ix(k);                                     % Current cell
       f = G.cells.faces(ixc(c) : ixc(c + 1) - 1, 1); % Faces on cell
       e = mcolon(ixf(f), ixf(f + 1) - 1);            % Edges on cell

       nodes  = unique(G.faces.nodes(e, 1));          % Unique nodes...
       coords = G.nodes.coords(nodes,:);            % ... and coordinates

       % Compute bounding box
       m = min(coords);
       M = max(coords);

       % Size of bounding box
       dx(k) = M(1) - m(1);
       if size(G.nodes.coords, 2) > 1,
          dy(k) = M(2) - m(2);
       else
          dy(k) = 1;
       end

       if size(G.nodes.coords, 2) > 2,
          dz(k) = G.cells.volumes(ix(k))/(dx(k)*dy(k));
       else
          dz(k) = 0;
       end
    end
end




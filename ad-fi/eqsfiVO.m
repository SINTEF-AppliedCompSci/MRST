function [eqs, state] = eqsfiVO(state0, state, dt, G, W, system, f, varargin)
% Generate equations for a Volatile 3Ph system (wet-gas, live-oil).

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

    opt = struct('Verbose',     mrstVerbose,...
                 'reverseMode', false,...
                 'scaling',     [],...
                 'resOnly',     false,...
                 'history',     [],  ...
                 'iteration',   -1,  ...
                 'stepOptions', []);

    opt = merge_options(opt, varargin{:});

    disgas = system.activeComponents.disgas;
    vapoil = system.activeComponents.vapoil;
    s = system.s;

    % current variables: ------------------------------------------------------
    p    = state.pressure;
    sW   = state.s(:,1);
    sG   = state.s(:,3);
    rs   = state.rs;
    rv   = state.rv;

    bhp = vertcat(state.wellSol.bhp);
    qWs    = vertcat(state.wellSol.qWs);
    qOs    = vertcat(state.wellSol.qOs);
    qGs    = vertcat(state.wellSol.qGs);

    % previous time-step variables ------------------------------------------------------
    p0  = state0.pressure;
    sW0 = state0.s(:,1);
    sG0 = state0.s(:,3);
    rs0 = state0.rs;
    rv0 = state0.rv;
    bhp0 = vertcat(state0.wellSol.bhp);
    qWs0 = vertcat(state0.wellSol.qWs);
    qOs0 = vertcat(state0.wellSol.qOs);
    qGs0 = vertcat(state0.wellSol.qGs);

    %Initialization of primary variables ----------------------------------
    [st1 , st2  , st3 ] = getCellStatus(state , disgas, vapoil);
    [st1p, st2p , st3p] = getCellStatus(state0, disgas, vapoil);
    if ~opt.resOnly,
        if ~opt.reverseMode,
            % define primary varible x and initialize
            x = st1.*rs + st2.*rv + st3.*sG;

            [p, sW, x, qWs, qOs, qGs, bhp] = ...
                initVariablesADI(p, sW, x, qWs, qOs, qGs, bhp);
            % define sG, rs and rv in terms of x
            sG = st2.*(1-sW) + st3.*x;
            if disgas
                rsSat = f.rsSat(p);
                rs = (~st1).*rsSat + st1.*x;
            else % otherwise rs = rsSat = const
                rsSat = rs;
            end
            if vapoil
                rvSat = f.rvSat(p);
                rv = (~st2).*rvSat + st2.*x;
            else % otherwise rv = rvSat = const
                rvSat = rv;
            end
        else
            x0 = st1p.*rs0 + st2p.*rv0 + st3p.*sG0;

            [p0, sW0, x0, zw, zw, zw, zw] = ...
                initVariablesADI(p0, sW0, x0, ...
                zeros(size(qWs0)) , zeros(size(qOs0)) , ...
                zeros(size(qGs0)) , zeros(size(bhp0)));                 %#ok
            sG0 = st2p.*(1-sW0) + st3p.*x0;
            if disgas
                rsSat0 = f.rsSat(p0);
                rs0 = (~st1p).*rsSat0  + st1p.*x0;
            else 
                rsSat0 = rs0; % Not used - remove
            end
            if vapoil
                rvSat0 = f.rvSat(p0);
                rv0 = (~st2p).*rvSat0  + st2p.*x0;
            else
                rvSat0 = rv0; % Not used - remove
            end
        end
    else % resOnly-case compute rsSat and rvSat for use in well eqs
        if disgas, rsSat = f.rsSat(p); else rsSat = rs; end
        if vapoil, rvSat = f.rvSat(p); else rvSat = rv; end
    end
    %----------------------------------------------------------------------
    %check for p-dependent tran mult:
    trMult = 1;
    if isfield(f, 'tranMultR'), trMult = f.tranMultR(p); end

    %check for p-dependent porv mult:
    pvMult = 1; pvMult0 = 1;
    if isfield(f, 'pvMultR')
        pvMult =  f.pvMultR(p);
        pvMult0 = f.pvMultR(p0);
    end

    %check for capillary pressure (p_cow)
    pcOW = 0;
    if isfield(f, 'pcOW')
        pcOW  = f.pcOW(sW);
    end
    %check for capillary pressure (p_cog)
    pcOG = 0;
    if isfield(f, 'pcOG')
        pcOG  = f.pcOG(sG);
    end

    % FLIUD PROPERTIES ---------------------------------------------------
    [krW, krO, krG] = f.relPerm(sW, sG);
    g  = norm(gravity);
    dz = s.grad(G.cells.centroids(:,3));

    % WATER PROPS (calculated at oil pressure)
    bW     = f.bW(p);
    rhoW   = bW.*f.rhoWS;
    % rhoW on face, avarge of neighboring cells (E100, not E300)
    rhoWf  = s.faceAvg(rhoW);
    mobW   = trMult.*krW./f.muW(p);
    dpW    = s.grad(p-pcOW) - g*(rhoWf.*dz);
    % water upstream-index
    upc  = (value(dpW)>=0);
    bWvW = s.faceUpstr(upc, bW.*mobW).*s.T.*dpW;

    % OIL PROPS
    if disgas
        bO  = f.bO(p, rs, ~st1);
        muO = f.muO(p, rs, ~st1);
    else
        bO  = f.bO(p);
        muO = f.muO(p);
    end
    rhoO   = bO.*(rs*f.rhoGS + f.rhoOS);
    rhoOf  = s.faceAvg(rhoO);
    mobO   = trMult.*krO./muO;
    dpO    = s.grad(p) - g*(rhoOf.*dz);
    % oil upstream-index
    upc = (value(dpO)>=0);
    bOvO   = s.faceUpstr(upc, bO.*mobO).*s.T.*dpO;
    if disgas, rsbOvO = s.faceUpstr(upc, rs).*bOvO;end

    % GAS PROPS (calculated at oil pressure)
    if vapoil
        bG  = f.bG(p, rv, ~st2);
        muG = f.muG(p, rv, ~st2);
    else
        bG  = f.bG(p);
        muG = f.muG(p);
    end
    rhoG   = bG.*(rv*f.rhoOS + f.rhoGS);
    rhoGf  = s.faceAvg(rhoG);
    mobG   = trMult.*krG./muG;
    dpG    = s.grad(p+pcOG) - g*(rhoGf.*dz);
    % gas upstream-index
    upc    = (value(dpG)>=0);
    bGvG   = s.faceUpstr(upc, bG.*mobG).*s.T.*dpG;
    if vapoil, rvbGvG = s.faceUpstr(upc, rv).*bGvG; end

    % EQUATIONS -----------------------------------------------------------
    sO  = 1- sW  - sG;
    sO0 = 1- sW0 - sG0;

    bW0 = f.bW(p0);
    if disgas, bO0 = f.bO(p0, rs0, ~st1p); else bO0 = f.bO(p0); end
    if vapoil, bG0 = f.bG(p0, rv0, ~st2p); else bG0 = f.bG(p0); end

    % oil eq:
    if vapoil
        eqs{1} = (s.pv/dt).*( pvMult.* (bO.* sO  + rv.* bG.* sG) - ...
                              pvMult0.*(bO0.*sO0 + rv0.*bG0.*sG0) ) + ...
                 s.div(bOvO + rvbGvG);
    else
        eqs{1} = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) + s.div(bOvO);
    end
    % water eq:
    eqs{2} = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + s.div(bWvW);
    % gas eq:
    if disgas
        eqs{3} = (s.pv/dt).*( pvMult.* (bG.* sG  + rs.* bO.* sO) - ...
                              pvMult0.*(bG0.*sG0 + rs0.*bO0.*sO0 ) ) + ...
                 s.div(bGvG + rsbOvO);
    else
        eqs{3} = (s.pv/dt).*( pvMult.*bG.*sG - pvMult0.*bG0.*sG0 ) + s.div(bGvG);
    end

    % well equations
    if ~isempty(W)
        wc    = vertcat(W.cells);
        if ~opt.reverseMode
            nperf = numel(wc);
            pw    = p(wc);
            rhows = [f.rhoWS, f.rhoOS, f.rhoGS];
            bw    = {bW(wc), bO(wc), bG(wc)};
            if ~disgas
                rsw = ones(nperf,1)*rs; rsSatw = ones(nperf,1)*rsSat; %constants
            else
                rsw = rs(wc); rsSatw = rsSat(wc);
            end
            if ~vapoil
                rvw = ones(nperf,1)*rv; rvSatw = ones(nperf,1)*rvSat; %constants
            else
                rvw = rv(wc); rvSatw = rvSat(wc);
            end
            rw    = {rsw, rvw};
            rSatw = {rsSatw, rvSatw};
            mw    = {mobW(wc), mobO(wc), mobG(wc)};

            optloc = {'iteration', opt.iteration, ...
                      'model', 'VO', ...
                      'allowWellSignChange', system.well.allowWellSignChange, ...
                      'allowControlSwitching', system.well.allowControlSwitching};
            
            [eqs(4:7), cqs, state.wellSol] = getWellContributions( W, state.wellSol, bhp, ...
                                                                      {qWs,qOs,qGs}, pw, rhows, ...
                                                                      bw, rw, rSatw, mw, ...
                                                                      optloc{:});

            eqs{1}(wc) = eqs{1}(wc) - cqs{2}; % Add src to oil eq
            eqs{2}(wc) = eqs{2}(wc) - cqs{1}; % Add src to water eq
            eqs{3}(wc) = eqs{3}(wc) - cqs{3}; % Add src to gas eq
        else
            % Force wells to be ADI variables.
            nw = numel(state.wellSol);
            zw = double2ADI(zeros(nw,1), p0);
            eqs(4:7) = {zw, zw, zw, zw};
        end
    else
        eqs(4:7) = {bhp, bhp, bhp, bhp};  % empty  ADIs
    end
end
%--------------------------------------------------------------------------
function [st1, st2, st3] = getCellStatus(state, disgas, vapoil)
% Status should be passed on from updateStateVO (to be sure definition is
% identical). rs and rv are assumed to be compatible, i.e. rx = rxSat for
% saturated cells and rx <= rxSat for undersaturated. Three values of
% status are:
% status 0: should not occur (almost water only -> state 3)
% status 1 oil, no gas  : x = rs, sg = 0    , rv = rvMax
% status 2 gas, no oil  : x = rv, sg = 1-sw , rs = rsMax
% status 3 oil and gas  : x = sg, rs = rsMax, rv = rvMax
if isfield(state, 'status')
    status = state.status;
else
    s = state.s;
    watOnly    = s(:,1) > 1- sqrt(eps);
    if ~vapoil
        oilPresent = true;
    else
        oilPresent = or(s(:,2) > 0, watOnly);
    end
    if ~disgas
        gasPresent = true;
    else
        gasPresent = or(s(:,3) > 0, watOnly);
    end
    status = oilPresent + 2*gasPresent;
end
if ~disgas
    st1 = false;
else
    st1 = status==1;
end
if ~vapoil
    st2 = false;
else
    st2 = status==2;
end
st3 = status == 3;
end


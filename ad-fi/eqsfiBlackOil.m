function [eqs, state, hst] = eqsfiBlackOil(state0, state, dt, G, W, system, f, varargin)
% Generate equations for a Black Oil system (oil, water and gas with gas dissolved in oil).

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


    if ~isempty(opt.scaling)
        scalFacs = opt.scaling;
    else
        scalFacs.rate = 1; scalFacs.pressure = 1;
    end

    hst = opt.history;
    s = system.s;

    % current variables: ------------------------------------------------------
    p    = state.pressure;
    sW   = state.s(:,1);
    sG   = state.s(:,3);
    rs   = state.rs;

    bhp = vertcat(state.wellSol.bhp);
    qWs    = vertcat(state.wellSol.qWs);
    qOs    = vertcat(state.wellSol.qOs);
    qGs    = vertcat(state.wellSol.qGs);

    % previous variables ------------------------------------------------------
    p0  = state0.pressure;
    sW0 = state0.s(:,1);
    sG0 = state0.s(:,3);
    rs0 = state0.rs;
    %--------------------------------------------------------------------------
    isSat = (sG>0) | (1 - sW - sG) < sqrt(eps);
    % Possible alternate condition...
    % isSat = (sG>0) | (1 - sW - sG) < sqrt(eps) | rs >= f.rsSat(p);

    if isempty(hst)
        hst.numch = zeros(G.cells.num,1);
    else
        hst.numch = hst.numch + double(xor(hst.isSat,isSat));
    end
    hst.isSat = isSat;

    %Initialization of independent variables ----------------------------------

    if ~opt.resOnly,
       % ADI variables needed since we are not only computing residuals.

       if ~opt.reverseMode,

          [p, sW, sG, rs, qWs, qOs, qGs, bhp] = ...
             initVariablesADI(p, sW, sG, rs, qWs, qOs, qGs, bhp);

       else

          [p0, sW0, sG0, rs0, zw, zw, zw, zw] = ...
             initVariablesADI(p0, sW0, sG0, rs0, ...
                              zeros(size(qWs)) , ...
                              zeros(size(qOs)) , ...
                              zeros(size(qGs)) , ...
                              zeros(size(bhp)));                      %#ok

       end
    end


    g  = norm(gravity);
    [Tw, dzw, Rw, wc, perf2well, pInx, iInxW] = getWellStuff(W);

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

    % -------------------------------------------------------------------------
    [krW, krO, krG] = f.relPerm(sW, sG);

    % water props (calculated at oil pressure OK?)
    bW     = f.bW(p);
    %bW     = f.bW(p-pcOW);
    rhoW   = bW.*f.rhoWS;
    % rhoW on face, avarge of neighboring cells (E100, not E300)
    rhoWf  = s.faceAvg(rhoW);
    mobW   = trMult.*krW./f.muW(p);
    %mobW   = trMult.*krW./f.muW(p-pcOW);
    dpW     = s.grad(p-pcOW) - g*(rhoWf.*s.grad(G.cells.centroids(:,3)));
    % water upstream-index
    upc = (double(dpW)>=0);
    bWvW = s.faceUpstr(upc, bW.*mobW).*s.T.*dpW;


    % oil props
    bO     = f.bO(p, rs, isSat);
    rhoO   = bO.*(rs*f.rhoGS + f.rhoOS);
    rhoOf  = s.faceAvg(rhoO);
    mobO   = trMult.*krO./f.muO(p,rs,isSat);
    dpO    = s.grad(p) - g*(rhoOf.*s.grad(G.cells.centroids(:,3)));
    % oil upstream-index
    upc = (double(dpO)>=0);
    bOvO   = s.faceUpstr(upc, bO.*mobO).*s.T.*dpO;
    rsbOvO = s.faceUpstr(upc, rs).*bOvO;

    % gas props (calculated at oil pressure OK?)
    bG     = f.bG(p);
    %bG     = f.bG(p+pcOG);
    rhoG   = bG.*f.rhoGS;
    rhoGf  = s.faceAvg(rhoG);
    mobG = trMult.*krG./f.muG(p);
    %mobG = trMult.*krG./f.muG(p+pcOG);

    dpG     = s.grad(p+pcOG) - g*(rhoGf.*s.grad(G.cells.centroids(:,3)));
    % water upstream-index
    upc = (double(dpG)>=0);
    bGvG = s.faceUpstr(upc,bG.*mobG).*s.T.*dpG;



    % EQUATIONS ---------------------------------------------------------------
    isSat0 = (double(sG0)>0);
    rsSat  = f.rsSat(p);


    % oil:
    eqs{1} = (s.pv/dt).*( pvMult.*bO.*(1-sW-sG) - pvMult0.*f.bO(p0,rs0,isSat0).*(1-sW0-sG0) ) + s.div(bOvO);
    % water:
    eqs{2} = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*f.bW(p0).*sW0 ) + s.div(bWvW);
    % gas:
    eqs{3} = (s.pv/dt).*...
        ( pvMult.*(bG.*sG + rs.*bO.*(1-sW-sG) ) -...
        pvMult0.*(f.bG(p0).*sG0 + rs0.*f.bO(p0,rs0,isSat0).*(1-sW0-sG0) ) )+ ...
        s.div(bGvG + rsbOvO);

    % closing eqs:
    eqs{4} = sG + 0*sG0;
    eqs{4}(isSat) = rs(isSat) - rsSat(isSat);

    % well equations
    if ~isempty(W)
        if ~opt.reverseMode
            pw   = p(wc);
            rhows = [f.rhoWS, f.rhoOS, f.rhoGS];
            %rhow = {rhoW(wc), rhoO(wc), rhoG(wc)};
            bw    = {bW(wc), bO(wc), bG(wc)};
            rw    = {rs(wc)};
            rSatw = {rsSat(wc)};
            mw    = {mobW(wc), mobO(wc), mobG(wc)};
            %mixOs = 1-mixWs-mixGs;
            [eqs(5:8), cqs, state.wellSol] = getWellContributions(...
                W, state.wellSol, bhp, {qWs,qOs,qGs}, pw, rhows, bw, rw, rSatw, mw, ...
                'iteration', opt.iteration);

            eqs{1}(wc) = eqs{1}(wc) - cqs{2}; % Add src to oil eq
            eqs{2}(wc) = eqs{2}(wc) - cqs{1}; % Add src to water eq
            eqs{3}(wc) = eqs{3}(wc) - cqs{3}; % Add src to gas eq
        else
            % Force wells to be ADI variables.
            nw = numel(state0.wellSol);
            zw = double2ADI(zeros(nw,1), p0);
            eqs(5:8) = {zw, zw, zw, zw};
        end
    else
        eqs(5:8) = {bhp, bhp, bhp, bhp};  % empty  ADIs
    end
end



%--------------------------------------------------------------------------


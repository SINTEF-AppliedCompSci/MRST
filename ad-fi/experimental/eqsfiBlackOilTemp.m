function [eqs, state, hst] = eqsfiBlackOilTemp(state0, state, dt, G, W, s, f, varargin)
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
                 'stepOptions', []);

    opt = merge_options(opt, varargin{:});


    if ~isempty(opt.scaling)
        scalFacs = opt.scaling;
    else
        scalFacs.rate = 1; scalFacs.pressure = 1;
    end

    hst = opt.history;

    % current variables: ------------------------------------------------------
    p    = state.pressure;
    T  = state.T;
    sW   = state.s(:,1);
    sG   = state.s(:,3);
    rs   = state.rs;

    pBHP = vertcat(state.wellSol.bhp);
    qWs  = vertcat(state.wellSol.qWs);
    qOs  = vertcat(state.wellSol.qOs);

    qGs  = vertcat(state.wellSol.qGs);

    % previous variables ------------------------------------------------------
    p0  = state0.pressure;
    T0  = state0.T;
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

    zw = zeros(size(pBHP));
    if ~opt.resOnly,
       % ADI variables needed since we are not only computing residuals.

       if ~opt.reverseMode,

          [p, sW, sG, rs, T, qWs, qOs, qGs, pBHP] = ...
             initVariablesADI(p, sW, sG, rs, T, qWs, qOs, qGs, pBHP);

       else

          [p0, sW0, sG0, rs0, zw, zw, zw, zw] = ...
             initVariablesADI(p0, sW0, sG0, rs0, T0,...
                              zeros(size(qGs)) , ...
                              zeros(size(qWs)) , ...
                              zeros(size(qOs)) , ...
                              zeros(size(pBHP)));                      %#ok

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
    mobW   = trMult.*krW./f.muW(p,T);
    %mobW   = trMult.*krW./f.muW(p-pcOW);
    dpW     = s.grad(p-pcOW) - g*(rhoWf.*s.grad(G.cells.centroids(:,3)));
    % water upstream-index
    upc = (double(dpW)>=0);
    bWvW = s.faceUpstr(upc, bW.*mobW).*s.T.*dpW;


    % oil props
    bO     = f.bO(p, rs, isSat);
    rhoO   = bO.*(rs*f.rhoGS + f.rhoOS);
    rhoOf  = s.faceAvg(rhoO);
    mobO   = trMult.*krO./f.muO(p,rs,isSat,T);
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
    sO=1-sW-sG;
    sO0=1-sW0-sG0;

    bO0=f.bO(p0,rs0,isSat0);bW0=f.bW(p0);bG0=f.bG(p0);
    % oil:
    eqs{1} = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) + s.div(bOvO);
    % water:
    eqs{2} = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + s.div(bWvW);
    % gas:
    eqs{3} = (s.pv/dt).*...
             ( pvMult.*(bG.*sG + rs.*bO.*sO ) -...
               pvMult0.*(bG0.*sG0 + rs0.*f.bO(p0,rs0,isSat0).*sO0 ) )+ ...
             s.div(bGvG + rsbOvO);

    % closing eqs:
    eqs{4} = sG + 0*sG0;
    eqs{4}(isSat) = rs(isSat) - rsSat(isSat);

    % add energy equations
    rhoS={f.rhoOS,f.rhoWS, f.rhoGS};
    bFsF={bO.*sO,bW.*sW,bG.*sG+rs.*bO.*sO};
    bFsF0={bO0.*sO0,bW0.*sW0,bG0.*sG0+rs0.*f.bO(p0,rs0,isSat0).*sO0};
    %{
    rhoF={rhoO,rhoW,rhoG};
    %rhoF0={bO0.*(rs0*f.rhoGS + f.rhoOS),f.rhoWS.*bW0,f.rhoGS.*bG0};
    rhoF0={bO0.*f.rhoOS,f.rhoWS.*bW0,f.rhoGS.*bG0};
    bF={bO,bW,bG}; bF0={bO0,bW0,bG0};
    sF={sO,sW,sG}; sF0={sO0,sW0,sG0};
    %}


    bFvF={bOvO,bWvW,bGvG+rsbOvO};

    eF={f.uO(p,T) , f.uW(p,T), f.uG(p,T)};
    eF0={f.uO(p0,T0) , f.uW(p0,T0), f.uG(p0,T0)};
    hF={f.hO(p,T), f.hW(p,T), f.hG(p,T)};

    uR=f.uR(T);uR0=f.uR(T0);
    vQ = s.T_r .* s.grad(T);
    wT=[W.T]';
    %vQqQ=W.WI_r*(wT-T(wc));
    vQqQ = 0*wT;
    %bFqF={bOqO,bWqW,bGqG};
    % well contributions is taken at the end

    
    vol=G.cells.volumes;
    eqs{5} = (1/dt).*((vol-pvMult.*s.pv).*uR-(vol-pvMult0.*s.pv).*uR0) + s.div( vQ);
    eqs{5}(wc) = eqs{5}(wc)-vQqQ;
    for i=1:numel(eF)       
        eqs{5}  =  eqs{5} + ((s.pv/dt).*( pvMult.*eF{i}.*rhoS{i}.*bFsF{i} - pvMult0.*eF0{i}.*rhoS{i}.*bFsF0{i} )...
                +  s.div( s.faceUpstr(bFvF{i}>0, rhoS{i}.*hF{i}) .* bFvF{i}));         
    end


    %solve local well-eqs
    if opt.stepOptions.solveWellEqs && ~opt.reverseMode
        state.wellSol = solveLocalWellEqs(W, pBHP, {qWs, qOs, qGs}, p(wc), ...
                                          {rhoW(wc), rhoO(wc), rhoG(wc)}, ...
                                          {bW(wc), bO(wc), bG(wc)}, rs(wc), ...
                                          {mobW(wc), mobO(wc), mobG(wc)}, ...
                                          state.wellSol, opt.stepOptions.maxitWell);
        if isa(pBHP, 'ADI')
            pBHP.val = vertcat(state.wellSol.bhp);
            qWs.val  = vertcat(state.wellSol.qWs);
            qOs.val  = vertcat(state.wellSol.qOs);
            qGs.val  = vertcat(state.wellSol.qGs);
        else
            pBHP = vertcat(state.wellSol.bhp);
            qWs  = vertcat(state.wellSol.qWs);
            qOs  = vertcat(state.wellSol.qOs);
            qGs  = vertcat(state.wellSol.qGs);
        end
    end

    % well equations/contributions
    [eqs(6:9), bq, state.wellSol] = ...
        getWellContributionsBO(W, pBHP, {qWs, qOs, qGs}, p(wc), ...
                                   {rhoW(wc), rhoO(wc), rhoG(wc)}, ...
                                   {bW(wc), bO(wc), bG(wc)}, rs(wc), ...
                                   {mobW(wc), mobO(wc), mobG(wc)}, ...
                                   state.wellSol);

    if opt.reverseMode
       % Force wells to be ADI variables.
       zeroW = 0*zw;
       for i = 5:8
          eqs{i} = eqs{i} + zeroW;
       end
    end

    eqs{1}(wc) = eqs{1}(wc) - bq{2};                 %oil
    eqs{2}(wc) = eqs{2}(wc) - bq{1};                 %water
    eqs{3}(wc) = eqs{3}(wc) - bq{3} - rs(wc).*bq{2}; %gas

    % wells for temprature
    % lett all perforation have entalpy
    bFqF={bq{2},bq{1},bq{3}};
    hFwp=cell(3,1);
    for i=1:numel(W)
       hFwp{1}=[hFwp{1};f.rhoOS*repmat(W(i).hO,numel(W(i).cells),1)];
       hFwp{2}=[hFwp{2};f.rhoWS*repmat(W(i).hW,numel(W(i).cells),1)];
       hFwp{3}=[hFwp{3};f.rhoGS*repmat(W(i).hG,numel(W(i).cells),1)];
    end
    HqH=cell(3,1);
    for i=1:3
        HqH{i}  = rhoS{i}.*hF{i}(wc).*bFqF{i};
        ind=bFqF{i}>0;
        HqH{i}(ind)  = hFwp{i}(ind).*bFqF{i}(ind);
    end
    %% add well contributioin
    for i=1:3
       eqs{5}(wc) = eqs{5}(wc)  -  HqH{i};
    end
    eqs{5}=eqs{5}/1e6;

end

%--------------------------------------------------------------------------


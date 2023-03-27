function [eqs, hst, explTrms] = eqsfiBlackOilExplicitWellsOGVE(state0, state, dt, G, W, s, f, varargin)
% Generate equations for a simplified Black Oil system (water and gas with gas dissolved in water).
opt = struct('Verbose',     mrstVerbose,...
             'reverseMode', false,...
             'scaling',     [],...
             'resOnly',     false,...
             'history',     [],...
             'bc',[]);
opt = merge_options(opt, varargin{:});

if ~isempty(opt.scaling)
    scalFacs = opt.scaling;
else
    scalFacs.rate = 1; scalFacs.pressure = 1;
end

hst = opt.history;

% current variables: ------------------------------------------------------
p     = state.pressure;
sG    = state.s(:,2);
sGmax = state.sGmax;
rs    = state.rs;
pBHP  = vertcat(state.wellSol.bhp);
qWs   = vertcat(state.wellSol.qWs);
qGs   = vertcat(state.wellSol.qGs);

% previous variables ------------------------------------------------------
p0     = state0.pressure;
sG0    = state0.s(:,2);
sGmax0 = state0.sGmax;

assert(all(state0.s(:)>=0))
assert(all(state.s(:)>=0))
rs0 = state0.rs;
%--------------------------------------------------------------------------
rsSat  = f.rsSat(p);
if(isfield(f,'dis_rate'))
    isSat = (sG>0) | rs>rsSat;% | (1 - sW + sG)  == 0;
else
    isSat = (sG>sqrt(eps)) | rs>rsSat;% | (1 - sW + sG)  == 0;
end
if isempty(hst)
    hst.numch = zeros(G.cells.num,1);
else
    hst.numch = hst.numch + double(xor(hst.isSat,isSat));
end
hst.isSat = isSat;

%Initialization of independent variables ----------------------------------
zw = zeros(size(pBHP));
if opt.resOnly
    % ADI variables aren't needed since we are only computing the residual.
    
elseif ~opt.reverseMode
    [p, sG, rs, qWs, qGs, pBHP,sGmax]  = initVariablesADI(p, sG, rs, qWs, qGs, pBHP,sGmax);
else
    [p0, sG0, rs0, ~, ~, zw] = initVariablesADI(p0, sG0, rs0,...
        zeros(size(qGs)),...
        zeros(size(qWs)), ...
        zeros(size(pBHP)));
end


g  = norm(gravity);
if(~isempty(W))
    [Tw, dzw, Rw, wc, perf2well, pInx, iInxG, iInxW] = getWellStuffWG(W);
end
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

%check for capillary pressure (p_cog)
pcWG = 0;
if isfield(f, 'pcWG')
    pcWG  = f.pcWG(sG,p,'sGmax',sGmax);
    pG=p;
    pG0=p0;
end

% -------------------------------------------------------------------------

krW = f.krW(1-sG, p,'sGmax',sGmax);
krG = f.krG(sG, p,'sGmax',sGmax);
% water props (calculated at oil pressure OK?)

% oil props
bW    = f.bW(p, rs, isSat);
rhoW  = bW.*(rs*f.rhoGS + f.rhoWS);
rhoWf = s.faceAvg(rhoW);
muW   = f.muW(p,rs,isSat);
mobW  = trMult.*krW./muW;

if(~any(strcmp(G.type,'topSurfaceGrid')))
    dpW    = s.grad(p) - g*(rhoWf.*s.grad(G.cells.centroids(:,3)));
    
else
    dpW    = s.grad(p) - g*(rhoWf.*s.grad(G.cells.z));
end

% oil upstream-index
upc = (double(dpW)>=0);
bWvW   = s.faceUpstr(upc, bW.*mobW).*s.T.*dpW;
rsbWvW = s.faceUpstr(upc, rs).*bWvW;

% gas props (calculated at oil pressure OK?)
bG    = f.bG(pG);
rhoG  = bG.*f.rhoGS;
rhoGf = s.faceAvg(rhoG);
muG   = f.muG(pG);%+pcWG);
mobG  = trMult.*krG./muG;

assert(all(double(rhoG)<double(rhoW)));

if(~any(strcmp(G.type,'topSurfaceGrid')))
    dpG     = s.grad(p+pcWG) - g*(rhoGf.*s.grad(G.cells.centroids(:,3)));
else
    dpG     = s.grad(p+pcWG) - g*(rhoGf.*s.grad(G.cells.z));
end

% gas upstream-index
upc = (double(dpG)>=0);
bGvG = s.faceUpstr(upc,bG.*mobG).*s.T.*dpG;

%WELLS ----------------------------------------------------------------
%mobWw  = mobW(wc);
%mobWw  = mobW(wc);
%mobGw  = mobG(wc);

% well bore head term computed explicitly, also compute connection phase
% volume fractions (alpha)

if(~isempty(W))
    pcWGw = pcWG(wc);
    if ~opt.resOnly
        rhow          = [rhoW.val(wc), rhoG.val(wc)];
    else
        rhow          = [rhoW(wc), rhoG(wc)];
    end
    [Hw, alpha]   = computeWellHeadOG(W, state.wellSol, rhow);
    
    % pressure drawdown
    
    %isInj    = double(drawdown) < 0;
    
    % connection mobilities
    bWw = bW(wc);
    bGw = bG(wc);
    rsw = rs(wc);
    pw  = p(wc);

    %set water injector mobility: mobw = mobw+mobo+mobg, mobo = 0;
    mobGw  = mobG(wc);
    mobWw  = mobW(wc);
    Tdrawdown = -Tw.*(-pBHP(perf2well)  + pw);
    %assert(all(~(iInxG& iInxW)));
    assert(isempty(intersect(iInxG, iInxW)));
    
    mobGw(iInxG) = mobG(iInxG)+mobW(iInxG);
    mobWw(iInxG) = 0;
    mobGw(iInxW) = 0;
    mobWw(iInxW) = mobG(iInxW)+mobW(iInxW);
    % assume all influx from a production well is water/(in represented as oil phase)
    pInxIn = pInx(Tdrawdown(pInx) < 0);
    mobGw(pInxIn) = 0;
    mobWw(pInxIn) = mobG(pInxIn)+mobW(pInxIn);
    rsw(pInxIn) = 0;
    
    assert(all(mobWw>=0))
    assert(all(mobGw>=0))
    
    qW = mobWw.*Tdrawdown;
    qG = mobGw.*Tdrawdown;
    %{
    qW = (-Tw).*mobWw.*(-(pBHP(perf2well) + 0.0*Hw) + p(wc));
    qG = (-Tw).*mobGw.*(-(pBHP(perf2well) + 0.0*Hw+0.0*pcWGw) + p(wc));
    %}
    bWqW  = bWw.*qW;
    bGqG  = bGw.*qG;
    
    
    % Compute explicit terms
    explTrms.wellFlux = [double(qW), double(qG)];
    
    if ~all(sign(Rw'*sum(explTrms.wellFlux,2))==vertcat(W.sign))
    %         warning('Some producers are becoming injectors or vice versa !!!!!!')
    %         fprintf('Wells changing roles...\n')
    end
else
    explTrms.wellFlux = [];%[double(qW), double(qG)];
end

if(~isempty(opt.bc))
    assert(all(strcmp(opt.bc.type,'pressure')));
    Tbc=s.T_all(opt.bc.face);
    bc_cell=sum(G.faces.neighbors(opt.bc.face,:),2);

    % assume cappillary pressure zero outside
    dzbc=(G.cells.z(bc_cell)-G.faces.z(opt.bc.face));
    pWbc=p(bc_cell);rhoWbc=rhoW(bc_cell);
    dpbc_o=opt.bc.value-pWbc+g*(rhoWbc.*dzbc);
    pGbc= pWbc+pcWG(bc_cell);rhoGbc=rhoG(bc_cell);
    dpbc_g=opt.bc.value-(pGbc)+g*(rhoGbc.*dzbc);
    bGmobGbc = bG(bc_cell).*mobG(bc_cell);
    bWmobWbc = bW(bc_cell).*mobW(bc_cell);
    bGmobGbc(dpbc_g>0)=0;
    rsbc    = rs(bc_cell);

    if(any(dpbc_o>0))
        bWmobWbc(dpbc_o>0)=bW(bc_cell(dpbc_o>0)).*(mobW(bc_cell(dpbc_o>0))+mobG(bc_cell(dpbc_o>0)));
        rsbc(dpbc_o>0) = 0;
    end
    bGqGbc  = -bGmobGbc.*Tbc.*(dpbc_g);
    bWqWbc  = -bWmobWbc.*Tbc.*(dpbc_o);
else
    bc_cell=[];
end


% EQUATIONS ---------------------------------------------------------------
isSat0 = (double(sG0)>0);
rsSat  = f.rsSat(p);


% oil:
eqs{1} = (s.pv/dt).*( pvMult.*bW.*(1-sG) - pvMult0.*f.bW(p0,rs0,isSat0).*(1-sG0) ) + s.div(bWvW);
if(~isempty(W))
    %[wc, cqs] = checkForRepititions(wc, {bWqW,bGqG});
    %bWqW=cqs{1};
    %bGqG=cqs{2};
    eqs{1}(wc) = eqs{1}(wc) - bWqW;
    %eqs{1}(wc(iInxW)) = eqs{1}(wc(iInxW)) - bWqW(iInxW);
    %eqs{1}(wc(pInx)) = eqs{1}(wc(pInx)) - bWqW(pInx);
end
% gas:
eqs{2} = (s.pv/dt).*...
    ( pvMult.*(bG.*sG + rs.*bW.*(1-sG) ) -...
    pvMult0.*(f.bG(pG0).*sG0 + rs0.*f.bW(p0,rs0,isSat0).*(1-sG0) ) )+ ...
    s.div(bGvG + rsbWvW); 

if(~isempty(W))
    %%{
    eqs{2}(wc(iInxG)) = eqs{2}(wc(iInxG)) - bGqG(iInxG);
    eqs{2}(wc(pInx)) = eqs{2}(wc(pInx)) - bGqG(pInx);%- rsw(pInx).*bWqW(pInx);
    pInxW=pInx(bWqW(pInx) > 0);% produce may have influx of water
    eqs{2}(wc(pInxW)) = eqs{2}(wc(pInxW))- rsw(pInxW).*bWqW(pInxW);
    %pInxW=pInx & bWqW > 0;
    %}
    %eqs{2}(wc) = eqs{2}(wc) - bGqG(wc) rsw(pInxW).*bWqW(pInxW);
end
if(~isempty(opt.bc))
    eqs{1}(bc_cell)  = eqs{1}(bc_cell) + bWqWbc;
    eqs{2}(bc_cell)  = eqs{2}(bc_cell) + bGqGbc+ rsbc.*bWqWbc;%could add dissolved gas
    assert(all(double(bGqGbc+ rsbc.*bWqWbc))>=0)
end

% closing eqs:
if(isfield(f,'dis_rate'))

    dis_rate=f.dis_rate.*(s.pv./G.cells.H);
    
    % set rate to zero if already saturated, or if there is no gas phase present
    
    % dis_rate adjustment, approach 1
    %meps=sqrt(eps);
    meps=eps*1000;
    dis_rate=dis_rate.*double((double(rs)<=(double(rsSat)-meps)) & (double(sG)>meps));
    %%{
    a=200;
    tanhyp=@(x,a) ((exp(a*x)-exp(-a*x))./(exp(a*x)+exp(-a*x)));
    s_fac=tanhyp(sG,a); % approximately one, but goes to 0 for very small values of sG
    rs_eps=(rsSat-rs)./f.dis_max;
    rs_fac=tanhyp(rs_eps,a);
    dis_rate=dis_rate.*s_fac.*rs_fac; % smoothly turn down dissolution rate when
                                      % sG goes to 0, or when dissolved value
                                      % approaches maximum.
    %}
    eqs{3} = (s.pv/dt).*...
        ( pvMult.*(rs.*bW.*(1-sG) ) -...
        pvMult0.*(rs0.*f.bW(p0,rs0,isSat0).*(1-sG0) ) )+ ...
        s.div(rsbWvW); 

    if(~isempty(W))
        eqs{3}(wc(pInx)) = eqs{3}(wc(pInx)) - rsw(pInx).*bWqW(pInx);
    end
    if(~isempty(opt.bc))
        eqs{3}(bc_cell)  = eqs{3}(bc_cell) + rsbc.*bWqWbc; 
    end
    eqs{3}=eqs{3}-dis_rate;
    
    
    %% Ensure dissolved quantity doesn't go below the minimum allowed
    % (i.e. that suggested by dissolution in zones with residual saturation)
    min_rs  = minRs(p,sG,sGmax,f,G);
    min_rs0 = minRs(p0,sG0,sGmax0,f,G);

    tmp = (s.pv/dt).*...
    ( double(pvMult).*(double(min_rs).*double(bW) ) -...
    pvMult0.*(rs0.*f.bW(p0,rs0,isSat0).*(1-sG0) ) )+ ...
    s.div(double(rsbWvW)); 
    
    if (~isempty(W))
        tmp(wc(pInx)) = tmp(wc(pInx)) - double(rsw(pInx)).*double(bWqW(pInx));
    end
    if(~isempty(opt.bc))
        tmp(bc_cell)  = tmp(bc_cell) + double(rsbc).*double(bWqWbc); 
    end
    tmp=tmp-double(dis_rate);
    ind_low_rs = tmp>-sqrt(eps);  % If so, then min_rs is larger than the
                                  % solution of eqs{3}, in other words, the
                                  % solution of eqs{3} is smaller than the
                                  % allowed value.  We have to modify eqs{3}.

    ind = ind_low_rs;  % Identify cells where 'rs' is moving below the minimum
                       % allowed value
    if(any(ind))
        % Force value of 'rs' in these cell to equal 'min_rs'
        eqs{3}(ind)=rs(ind).*(1-sG(ind))-min_rs(ind); % @@ multiply by bW?????
        eqs{3}(ind)=eqs{3}(ind).*(s.pv(ind)/dt);
    end
    
    %% Ensure dissolved quantity doesn't exceed 'rsSat'
    
    is_sat_loc=((double(rs)>=double(rsSat)) &  (eqs{3}<0)) ;
    % force saturation in these cells to equal 'rsSat'.
    eqs{3}(is_sat_loc) = (rs(is_sat_loc) - rsSat(is_sat_loc)).*s.pv(is_sat_loc)/dt;

    %% Compute changes to sGmax

    % Default equation: force 'sGmax' to equal 'sG'
    eqs{7}=(sGmax-sG).*(s.pv/dt);
    tmp=(s.pv/dt).*...
        (pvMult.*bG.*sGmax- pvMult.*f.bG(pG).*sGmax0)*f.res_gas./(1-f.res_water)+dis_rate;
    tmp2=(s.pv/dt).*...
         (pvMult.*double(bG).*double(sG)- pvMult.*f.bG(pG).*sGmax0)*f.res_gas./(1-f.res_water)+dis_rate;
    % tmp=(s.pv/dt).*...
    %     (pvMult.*bG.*sGmax- pvMult0.*f.bG(pG0).*sGmax0)*f.res_gas./(1-f.res_water)+dis_rate;
    % tmp2=(s.pv/dt).*...
    %      (pvMult.*double(bG).*double(sG)- pvMult0.*f.bG(pG0).*sGmax0)*f.res_gas./(1-f.res_water)+dis_rate;
    %%{
    % the case of disolution do not manage to remove all residual
    % saturation  and the new state do not is not fully dissolved
    %ind1= (tmp2<0) & (sGmax < sGmax0) & (~is_sat_loc);
    ind1= (tmp2<0) & (~is_sat_loc);
    if(any(ind1))
        eqs{7}(ind1)=tmp(ind1);
    end
    % the new state is fully dissolved but all residual saturation is
    % not gone
    ind3= is_sat_loc & (sGmax < sGmax0) & (sGmax > sG);
    if(any(ind3))
        tmp = (s.pv/dt).*...
        (pvMult.*bG.*sGmax- pvMult.*f.bG(pG).*sGmax0)*f.res_gas./(1-f.res_water)+(rs-rs0);
        %(pvMult.*bG.*sGmax- pvMult0.*f.bG(pG0).*sGmax0)*f.res_gas./(1-f.res_water)+(rs-rs0);
        eqs{7}(ind3)= tmp(ind3);
    end
    
else
    eqs{3} = sG + 0*sG0;
    eqs{3}(isSat) = rs(isSat) - rsSat(isSat);
    eqs{7}=sGmax-max(sGmax0,sG);
end

% Well equations
% Force wells to be independent ADI variables


% Last eq: boundary cond
zeroW=0*pBHP;
if(~isempty(W))
    eqs{4} = -Rw'*bWqW + qWs + zeroW;
    eqs{5} = -Rw'*(bGqG + rsw.*bWqW) + qGs + zeroW;
    eqs{6} = handleBC(W, pBHP, [], qWs, qGs, scalFacs) + zeroW;
else
    eqs{4} = zeroW;
    eqs{5}  = zeroW;
    eqs{6} = zeroW;
end
for i=[1,2,3,7]
    eqs{i}=eqs{i}*dt/year;
end

end

%--------------------------------------------------------------------------

function [Hw, alpha] = computeWellHeadOG(W, wellSol, rho)
inx = 0;
if ~isfield(wellSol, 'flux') || isempty(wellSol(1).flux)
    for wnr = 1:numel(W)
        wellSol(wnr).flux = ones(numel(W(wnr).cells), 1)*W(wnr).compi(:,[2 3])*W(wnr).sign;
    end
end
rhoMix = cell(numel(W),1);
alpha  = cell(numel(W),1);

% the following could be replaced by e.g. a linear system for multilateral
% wells
for wnr = 1:numel(W)
    wbOut = [0 0];
    alpha{wnr}  = zeros(numel(W(wnr).cells), 2);
    rhoMix{wnr} = zeros(numel(W(wnr).cells), 1);
    for cnr = numel(W(wnr).cells):-1:1
        fluxOut = wellSol(wnr).flux(cnr,:);
        if cnr == 1 && W(wnr).sign==1
            wbIn    = sum(wbOut + fluxOut)*W(wnr).compi(:,[2 3]);
        else
            wbIn    = wbOut + fluxOut;
        end
        totIn   = (wbIn    > 0).*wbIn  + ...
            (wbOut   < 0).*wbOut + ...
            (fluxOut < 0).*fluxOut;
        if abs(sum(totIn)) < 1e-6/day
            totIn = W(wnr).compi(:,[2 3]);
        end
        if(sum(totIn)>0)
            a = totIn/sum(totIn);
        else
            a=W(wnr).compi(:,[2 3]);
        end
        alpha{wnr}(cnr,:)  = a;
        rhoMix{wnr}(cnr) = sum(a.*rho(inx+cnr,:));
        
        wbOut = wbIn;
    end
    inx = inx + numel(W(wnr).cells);
end
rhoMixAvg = cellfun(@(x)[x(1); .5*(x(1:end-1)+x(2:end))], rhoMix, 'UniformOutput', false);

g   = norm(gravity);
dzw = arrayfun(@(x)[x.dZ(1); x.dZ(2:end)-x.dZ(1:end-1)], W, 'UniformOutput', false);
ddp = cellfun(@(x,y)g*x.*y, dzw, rhoMixAvg, 'UniformOutput', false);
Hw  = cellfun(@cumsum, ddp, 'UniformOutput', false);
Hw = vertcat(Hw{:});
alpha = vertcat(alpha{:});
%dzw = vertcat(W.dZ);
%dpw = g*dzw.*vertcat(rhoMixAvg{:});
end

%--------------------------------------------------------------------------

function [mO, mG, crossFlow] = computeWellMobilitiesOG(W, qWs, qGs, mO,  mG,  bW, bG, rs, isInj)
% for producer producing connections, m = [mW mO mG] remains unchanged
% for producer injecting connections the (assumed uniform) flow in the
% well-bore is calculated by
%   qW = qWs/bW
%   qW = qWs/bW
%   qG = (qGs-rs*qWs)/bG
% and the mobility as m = [qW qW qG]*(mW+mO+mG)/(qW+qW+qG)
%
% for injector injecting connections the mobility is calculated by
% m = compi*(mW+mO+mG)
% injector producing connections is treated as above, but supresses a
% warning

nPerf = arrayfun(@(x)numel(x.cells), W)';
inj   = vertcat(W.sign)==1;
iInx  = rldecode(inj, nPerf);
% update crossflow flag
%check for injector producing connections
% if any(and(iInx, ~isInj))
%     warning('Crossflow detected in injectors, no special treatment for this case')
% end
if(any(iInx))
    injComp = rldecode(vertcat(W(inj).compi), nPerf(inj));
    mtInj   = mO(iInx) + mG(iInx);
    mO(iInx) = injComp(:,2).*mtInj;
    mG(iInx) = injComp(:,3).*mtInj;
end

crossFlow = and(~iInx,isInj);
if any(crossFlow)
    ix = find(crossFlow);
    perf2well = rldecode((1:numel(W))', nPerf);
    qW = qWs(perf2well(ix))./bW(ix);
    qG = (qGs(perf2well(ix))-rs(ix).*qWs(perf2well(ix)))./bG(ix);
    qt =  qW + qG;
    
    % don't bother with connections having almost zero flow
    %     zeroInx = abs(qt.val)<1e-6/day;
    %     ix = ix(~zeroInx);
    if ~isempty(ix)
        fprintf('Crossflow detected in %2.0d connections in wells ', numel(ix));
        fprintf('%d ', unique(perf2well(ix)))
        fprintf('\n')
        mt =  mO(ix) + mG(ix);
        mO(ix) = qW.*mt./qt;
        mG(ix) = qG.*mt./qt;
        
    end
end
assert(all(double(mO)>=0))
assert(all(double(mG)>=0))
end
function [wc, cqs] = checkForRepititions(wc, cqs)
[c, ia, ic] = unique(wc);%, 'stable');
if numel(c) ~= numel(wc)
    A = sparse(ic, (1:numel(wc))', 1, numel(c), numel(wc));
    wc = c;
    for k=1:numel(cqs)
        cqs{k} = A*cqs{k};
    end
end
end



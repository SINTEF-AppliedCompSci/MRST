function [eqs, hst] = eqsfiOGExplicitWellsVE(state0, state, dt, G, W, system, f, varargin)

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'scaling', [],...
             'resOnly', false,...
             'history', [],...
             'bc',[],...
             'iteration',-1,...% not used but defined to avoid warning in adjoint simulations
             'stepOptions',[]); %compatibility only

opt = merge_options(opt, varargin{:});

if ~isempty(opt.scaling)
    scalFacs = opt.scaling;
else
    scalFacs.rate = 1; scalFacs.pressure = 1;
end
s=system;
hst = opt.history;

% current variables: ------------------------------------------------------
p    = state.pressure;
sG   = state.s(:,2);
if(isfield(state0,'smax'))
    sGmax=max(state.s(:,2),state0.smax(:,2));
    sGmax0=state0.smax(:,2);
else
    sGmax=[];
end
pBHP = vertcat(state.wellSol.bhp);
qGs  = vertcat(state.wellSol.qGs);
qWs  = vertcat(state.wellSol.qWs);

% previous variables ------------------------------------------------------
p0  = state0.pressure;
sG0 = state0.s(:,2);
%--------------------------------------------------------------------------


%Initialization of independent variables ----------------------------------

zw = zeros(size(pBHP));
if opt.resOnly
    % ADI variables aren't needed since we are only computing the residual.
elseif ~opt.reverseMode
    [p, sG, qGs, qWs, pBHP]  = initVariablesADI(p, sG, qGs, qWs, pBHP);
else
    [p0, sG0, ~, ~, zw] = initVariablesADI(p0, sG0,...
                                            zeros(size(qGs)), ...
                                            zeros(size(qWs)), ...
                                            zeros(size(pBHP))...
                                            );
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



% -------------------------------------------------------------------------
% [krW, krW] = f.relPerm(sW);
krG = f.krG(sG, p, 'sGmax',sGmax);
krW = f.krW(1-sG, p, 'sGmax',sGmax);
%check for capillary pressure (p_cow)
pcWG = 0;
if isfield(f, 'pcWG') 
    pcWG  = f.pcWG(sG, p,'sGmax',sGmax);
    %pcWG0  = f.pcWG(sG0, p0,'sGmax',sGmax0);
    % this is used for all density evaluations
    pG=p;%+pcWG;
    pG0=p0;%+pcWG0;%+
    
else
   pG=p;%+pcWG;
   pG0=p0; 
end
% gas props 

bG     = f.bG(pG);
rhoG   = bG.*f.rhoGS;

rhoGf  = s.faceAvg(rhoG);

mobG   = trMult.*krG./f.muG(pG);
if(~any(strcmp(G.type,'topSurfaceGrid')))
    dpG     = s.grad(p+pcWG) - g*(rhoGf.*s.grad(G.cells.centroids(:,3)));
else
    dpG     = s.grad(pG+pcWG) - g*(rhoGf.*s.grad(G.cells.z));
end
% gas upstream-index
upc = (double(dpG)>=0);
bGvG = s.faceUpstr(upc, bG.*mobG).*s.T.*dpG;


% water props
bW     = f.bW(p);
rhoW   = bW.*f.rhoWS;
rhoWf  = s.faceAvg(rhoW);
%mobW   = trMult.*krW./f.BWxmuW(p);
mobW   = trMult.*krW./(f.BW(p).*f.muW(p));
if(~any(strcmp(G.type,'topSurfaceGrid')))
    dpW    = s.grad(p) - g*(rhoWf.*s.grad(G.cells.centroids(:,3)));
else
    dpW    = s.grad(p) - g*(rhoWf.*s.grad(G.cells.z));
end
% oil upstream-index
upc = (double(dpW)>=0);
bWvW   = s.faceUpstr(upc, bW.*mobW).*s.T.*dpW;


%WELLS ----------------------------------------------------------------
if(~isempty(W))
    pcWGw = pcWG(wc);
    bGw     = bG(wc);
    bWw     = bW(wc);
    mobGw  = mobG(wc);
    mobWw  = mobW(wc);

%producer mobility
bGmobGw  = bGw.*mobGw;
bWmobWw  = bWw.*mobWw;

%set water injector mobility: mobw = mobw+mobo+mobg, mobo = 0;
if(true)
    bGmobGw(iInxG) = bGw(iInxG).*(mobGw(iInxG) + mobWw(iInxG));
    bWmobWw(iInxG) = 0;
    bGmobGw(iInxW) = 0; 
    bWmobWw(iInxW) = bWw(iInxW).*(mobGw(iInxW) + mobWw(iInxW));
else
   bGmobGw(iInxG) = bGw(iInxG)./1e-3;
   bWmobWw(iInxG) = 0;
   bGmobGw(iInxW) = 0; 
   bWmobWw(iInxW) = bWw(iInxW)./1e-3;
end
   
pw  = p(wc);
perfP=pBHP(perf2well);
bGqG  = -bGmobGw.*Tw.*(perfP - pw + 0.0*pcWGw + 0.0*g*dzw.*rhoG(wc));
bWqW  = -bWmobWw.*Tw.*(perfP - pw + 0.0*g*dzw.*rhoW(wc));
end
% add contributions from boundary
if(~isempty(opt.bc))
    assert(all(strcmp(opt.bc.type,'pressure')));
    Tbc=s.T_all(opt.bc.face);
    bc_cell=sum(G.faces.neighbors(opt.bc.face,:),2);
    %bc_cell_loc=false(G.cells.num,1);
    %bc_cell_loc(bc_cell_nr)=true;
    % assume cappillary pressure zero outside
    dzbc=(G.cells.z(bc_cell)-G.faces.z(opt.bc.face));
    pWbc=p(bc_cell);rhoWbc=rhoW(bc_cell);    
    dpbc_o=opt.bc.value-pWbc+g*(rhoWbc.*dzbc);
    pGbc= pWbc+pcWG(bc_cell);rhoGbc=rhoG(bc_cell);
    dpbc_g=opt.bc.value-(pGbc)+g*(rhoGbc.*dzbc);
    bGmobGbc = bG(bc_cell).*mobG(bc_cell);
    bWmobWbc = bW(bc_cell).*mobW(bc_cell);
    bGmobGbc(dpbc_g>0)=0;
    %bGmobGbc(dpbc_o>0)=0;
    %bWmobWbc(dpbc_g>0)=0;
    if(any(dpbc_o>0))
        bWmobWbc(dpbc_o>0)=bW(bc_cell(dpbc_o>0)).*(mobW(bc_cell(dpbc_o>0))+mobG(bc_cell(dpbc_o>0)));
    end
    bGqGbc  = -bGmobGbc.*Tbc.*(dpbc_g);
    bWqWbc  = -bWmobWbc.*Tbc.*(dpbc_o);
else
    bc_cell=[];
end


% EQUATIONS ---------------------------------------------------------------


% oil:
eqs{1} = (s.pv/dt).*( pvMult.*bW.*(1-sG) - pvMult0.*f.bW(p0).*(1-sG0) ) + s.div(bWvW);


% water:
eqs{2} = (s.pv/dt).*( pvMult.*bG.*sG - pvMult0.*f.bG(pG0).*sG0 ) + s.div(bGvG);

    if(~isempty(W))
        [wc, cqs] = checkForRepititions(wc, {bWqW,bGqG});
        eqs{1}(wc) = eqs{1}(wc) + cqs{1};
        eqs{2}(wc) = eqs{2}(wc) + cqs{2};
        %eqs{1}(wc) = eqs{1}(wc) + bWqW;
        %eqs{2}(wc) = eqs{2}(wc) + bGqG;
    end
    if(~isempty(opt.bc))
        eqs{1}(bc_cell)  = eqs{1}(bc_cell) + bWqWbc;
        eqs{2}(bc_cell)  = eqs{2}(bc_cell) + bGqGbc;
        assert(sum(double(bGqGbc))>=0)
    end

if(~opt.reverseMode)
    % well equations
    zeroW = 0*pBHP;
    if(~isempty(W))
        eqs{3} = Rw'*bGqG + qGs + zeroW;
        eqs{4} = Rw'*bWqW + qWs + zeroW;
        
        % Last eq: boundary cond
        eqs{5} = handleBC(W, pBHP, [] , qWs, qGs, scalFacs) + zeroW;
    else
        eqs{3}=zeroW;
        eqs{4}=zeroW;
        eqs{5}=zeroW;
    end
    
else
    % in reverse mode just gather zero-eqs of correct size
    for eqn = 3:5
        nw = numel(state0.wellSol);
        zw = double2ADI(zeros(nw,1), p0);
        eqs(3:5) = {zw, zw, zw};
    end
end



for jj=1:numel(eqs)
   assert(all(isfinite(double(eqs{jj})))); 
end

%ol=double(ones(1,G.cells.num)*(s.pv.*f.pvMultR(p).*f.bG(p).*sG));% add mass later
%ouble(vol*f.rhoGS/1e9)
%sum([double(eqs{2})*dt,...
if(false)
if(~isempty(W))
    disp([sum(double(s.pv.*(pvMult.*bG.*sG)))-... % mass CO2 after
        sum(double(s.pv.*(pvMult0.*f.bG(pG0).*sG0))),...% mass CO2 before
        sum(double(s.div(bGvG)))*dt,...
        sum(double(bGqGbc))*dt,...
        sum(double(bGqG))*dt]*f.rhoG)% injected CO2
else
   disp([sum(double(s.pv.*(pvMult.*bG.*sG)))-... % mass CO2 after
        sum(double(s.pv.*(pvMult0.*f.bG(pG0).*sG0))),...% mass CO2 before
        sum(double(s.div(bGvG)))*dt,...
        sum(double(bGqGbc))*dt]);%,...
end
end
end
%--------------------------------------------------------------------------
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

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
s=system.s;
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
qOs  = vertcat(state.wellSol.qOs);

% previous variables ------------------------------------------------------
p0  = state0.pressure;
sG0 = state0.s(:,2);
%--------------------------------------------------------------------------


%Initialization of independent variables ----------------------------------

zw = zeros(size(pBHP));
if opt.resOnly
    % ADI variables aren't needed since we are only computing the residual.
elseif ~opt.reverseMode
    [p, sG, qGs, qOs, pBHP]  = initVariablesADI(p, sG, qGs, qOs, pBHP);
else
    [p0, sG0, ~, ~, zw] = initVariablesADI(p0, sG0,...
                                            zeros(size(qGs)), ...
                                            zeros(size(qOs)), ...
                                            zeros(size(pBHP))...
                                            );
end

g  = norm(gravity);
if(~isempty(W))
    [Tw, dzw, Rw, wc, perf2well, pInx, iInxG, iInxO] = getWellStuffOG(W);
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
% [krW, krO] = f.relPerm(sW);
krG = f.krG(sG, p, 'sGmax',sGmax);
krO = f.krOG(1-sG, p, 'sGmax',sGmax);
%check for capillary pressure (p_cow)
pcOG = 0;
if isfield(f, 'pcOG') 
    pcOG  = f.pcOG(sG, p,'sGmax',sGmax);
    %pcOG0  = f.pcOG(sG0, p0,'sGmax',sGmax0);
    % this is used for all density evaluations
    pG=p;%+pcOG;
    pG0=p0;%+pcOG0;%+
    
else
   pG=p;%+pcOG;
   pG0=p0; 
end
% water props (calculated at oil pressure OK?)
%bW     = f.bW(p);
bG     = f.bG(pG);
rhoG   = bG.*f.rhoGS;
% rhoW on face, avarge of neighboring cells (E100, not E300)
rhoGf  = s.faceAvg(rhoG);
%mobW   = trMult.*krW./f.muW(p);
mobG   = trMult.*krG./f.muG(pG);
if(~any(strcmp(G.type,'topSurfaceGrid')))
    dpG     = s.grad(p+pcOG) - g*(rhoGf.*s.grad(G.cells.centroids(:,3)));
else
    dpG     = s.grad(pG+pcOG) - g*(rhoGf.*s.grad(G.cells.z));
end
% water upstream-index
upc = (double(dpG)>=0);
bGvG = s.faceUpstr(upc, bG.*mobG).*s.T.*dpG;


% oil props
bO     = f.bO(p);
rhoO   = bO.*f.rhoOS;
rhoOf  = s.faceAvg(rhoO);
%mobO   = trMult.*krO./f.BOxmuO(p);
mobO   = trMult.*krO./(f.BO(p).*f.muO(p));
if(~any(strcmp(G.type,'topSurfaceGrid')))
    dpO    = s.grad(p) - g*(rhoOf.*s.grad(G.cells.centroids(:,3)));
else
    dpO    = s.grad(p) - g*(rhoOf.*s.grad(G.cells.z));
end
% oil upstream-index
upc = (double(dpO)>=0);
bOvO   = s.faceUpstr(upc, bO.*mobO).*s.T.*dpO;


%WELLS ----------------------------------------------------------------
if(~isempty(W))
    pcOGw = pcOG(wc);
    bGw     = bG(wc);
    bOw     = bO(wc);
    mobGw  = mobG(wc);
    mobOw  = mobO(wc);

%producer mobility
bGmobGw  = bGw.*mobGw;
bOmobOw  = bOw.*mobOw;

%set water injector mobility: mobw = mobw+mobo+mobg, mobo = 0;
if(true)
    bGmobGw(iInxG) = bGw(iInxG).*(mobGw(iInxG) + mobOw(iInxG));
    bOmobOw(iInxG) = 0;
    bGmobGw(iInxO) = 0; 
    bOmobOw(iInxO) = bOw(iInxO).*(mobGw(iInxO) + mobOw(iInxO));
else
   bGmobGw(iInxG) = bGw(iInxG)./1e-3;
   bOmobOw(iInxG) = 0;
   bGmobGw(iInxO) = 0; 
   bOmobOw(iInxO) = bOw(iInxO)./1e-3;
end
   
pw  = p(wc);
perfP=pBHP(perf2well);
bGqG  = -bGmobGw.*Tw.*(perfP - pw + 0.0*pcOGw + 0.0*g*dzw.*rhoG(wc));
bOqO  = -bOmobOw.*Tw.*(perfP - pw + 0.0*g*dzw.*rhoO(wc));
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
    pObc=p(bc_cell);rhoObc=rhoO(bc_cell);    
    dpbc_o=opt.bc.value-pObc+g*(rhoObc.*dzbc);
    pGbc= pObc+pcOG(bc_cell);rhoGbc=rhoG(bc_cell);
    dpbc_g=opt.bc.value-(pGbc)+g*(rhoGbc.*dzbc);
    bGmobGbc = bG(bc_cell).*mobG(bc_cell);
    bOmobObc = bO(bc_cell).*mobO(bc_cell);
    bGmobGbc(dpbc_g>0)=0;
    %bGmobGbc(dpbc_o>0)=0;
    %bOmobObc(dpbc_g>0)=0;
    if(any(dpbc_o>0))
        bOmobObc(dpbc_o>0)=bO(bc_cell(dpbc_o>0)).*(mobO(bc_cell(dpbc_o>0))+mobG(bc_cell(dpbc_o>0)));
    end
    bGqGbc  = -bGmobGbc.*Tbc.*(dpbc_g);
    bOqObc  = -bOmobObc.*Tbc.*(dpbc_o);
else
    bc_cell=[];
end


% EQUATIONS ---------------------------------------------------------------


% oil:
eqs{1} = (s.pv/dt).*( pvMult.*bO.*(1-sG) - pvMult0.*f.bO(p0).*(1-sG0) ) + s.div(bOvO);


% water:
eqs{2} = (s.pv/dt).*( pvMult.*bG.*sG - pvMult0.*f.bG(pG0).*sG0 ) + s.div(bGvG);

    if(~isempty(W))
        [wc, cqs] = checkForRepititions(wc, {bOqO,bGqG});
        eqs{1}(wc) = eqs{1}(wc) + cqs{1};
        eqs{2}(wc) = eqs{2}(wc) + cqs{2};
        %eqs{1}(wc) = eqs{1}(wc) + bOqO;
        %eqs{2}(wc) = eqs{2}(wc) + bGqG;
    end
    if(~isempty(opt.bc))
        eqs{1}(bc_cell)  = eqs{1}(bc_cell) + bOqObc;
        eqs{2}(bc_cell)  = eqs{2}(bc_cell) + bGqGbc;
        assert(sum(double(bGqGbc))>=0)
    end

if(~opt.reverseMode)
    % well equations
    zeroW = 0*pBHP;
    if(~isempty(W))
        eqs{3} = Rw'*bGqG + qGs + zeroW;
        eqs{4} = Rw'*bOqO + qOs + zeroW;
        
        % Last eq: boundary cond
        eqs{5} = handleBC(W, pBHP, [] , qOs, qGs, scalFacs) + zeroW;
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










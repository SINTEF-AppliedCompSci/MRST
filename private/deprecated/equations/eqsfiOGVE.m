function [eqs, hst] = eqsfiOGVE(state0, state, dt, G, W, system, f, varargin)

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
% water props (calculated at oil pressure OK?)
%bW     = f.bW(p);
bG     = f.bG(pG);
rhoG   = bG.*f.rhoGS;
% rhoW on face, avarge of neighboring cells (E100, not E300)
rhoGf  = s.faceAvg(rhoG);
%mobW   = trMult.*krW./f.muW(p);
mobG   = trMult.*krG./f.muG(pG);
if(~any(strcmp(G.type,'topSurfaceGrid')))
    dpG     = s.grad(p+pcWG) - g*(rhoGf.*s.grad(G.cells.centroids(:,3)));
else
    dpG     = s.grad(pG+pcWG) - g*(rhoGf.*s.grad(G.cells.z));
end
% water upstream-index
upc = (double(dpG)>=0);
bGvG = s.faceUpstr(upc, bG.*mobG).*s.T.*dpG;


% oil props
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



if ~isempty(W)
    if ~opt.reverseMode
        wc    = vertcat(W.cells);
        pw   = p(wc);
        rhos = [f.rhoWS, f.rhoGS];
        bw   = {bW(wc), bG(wc)};
        rw   = {};
        mw   = {mobW(wc), mobG(wc)};
        optloc = {'iteration', opt.iteration, ...
                  'model', 'OG', ...
                  'allowWellSignChange', system.well.allowWellSignChange, ...
                  'allowControlSwitching', system.well.allowControlSwitching};
        
        [eqs(3:5), cqs, state.wellSol] = getWellContributions(W, state.wellSol, pBHP, {qWs, qGs}, ...
                                                                 pw, rhos, bw, rw, rw, mw, ...
                                                                 optloc{:});

        [wc, cqs] = checkForRepititions(wc, cqs);
        eqs{1}(wc) = eqs{1}(wc) - cqs{1};
        eqs{2}(wc) = eqs{2}(wc) - cqs{2};
    else
        % in reverse mode just gather zero-eqs of correct size
        for eqn = 3:5
            nw = numel(state.wellSol);
            zw = double2ADI(zeros(nw,1), p0);
            eqs(3:5) = {zw, zw, zw};
        end
    end
else % no wells
    eqs(3:5) = {pBHP, pBHP, pBHP};  % empty  ADIs
end 



for jj=1:numel(eqs)
   assert(all(isfinite(double(eqs{jj})))); 
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










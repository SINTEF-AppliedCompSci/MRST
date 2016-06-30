function [problem, state] = equationsWaterThermal(state0, state, model, dt, drivingForces, varargin)
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

[p, T, wellSol] = model.getProps(state, 'pressure','temperature', 'wellsol');

[p0 , T0] = model.getProps(state0, 'pressure','temperature');

pBH    = vertcat(wellSol.bhp);
qWs    = vertcat(wellSol.qWs);

%Initialization of independent variables ----------------------------------

if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [p,T, qWs, pBH] = ...
            initVariablesADI(p, T, qWs, pBH);
    else
        [p0, T0, tmp, tmp] = ...
            initVariablesADI(p0, sW0,          ...
            zeros(size(qWs)), ...
            zeros(size(qOs)), ...
            zeros(size(pBH)));                          %#ok
    end
end


clear tmp
%grav  = gravity;
gdz   = s.Grad(G.cells.centroids) * model.gravity';
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
bW     = f.bW(p, T);bW0 = f.bW(p0,T0);
rhoW   = bW.*f.rhoWS;
% rhoW on face, avarge of neighboring cells (E100, not E300)
rhoWf  = s.faceAvg(rhoW);
mobW   = trMult./f.muW(p,T);
dpW     = s.Grad(p) - rhoWf.*gdz;
% water upstream-index
upcw = (double(dpW)<=0);
vW = - s.faceUpstr(upcw, mobW).*trans.*dpW;
bWvW = s.faceUpstr(upcw, bW).*vW;

% get themal propeties
uR = T.*model.rock.cR.*model.rock.rhoR;uR0 = T0.*model.rock.cR.*model.rock.rhoR;
uW = f.uW(p,T);uW0 = f.uW(p0,T0);
hW = f.hW(p,T);

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
eqs{1} = (s.pv/dt).*( pvMult.*bW - pvMult0.*bW0 ) + s.Div(bWvW);
[eqs, qBC, qRes, BCtocellMap, qSRC, srcCells, bcCells] = addFluxesFromSourcesAndBC(model, eqs, ...
                                          {p},...
                                        {rhoW},...
                                        {mobW}, ...
                                        {bW},  ...
                                       {ones(G.cells.num,1)}, ...
                                       drivingForces);
                                   
% accumulation of energy and conduction of heat in rock
vol=model.G.cells.volumes;
vQ = -s.T_r .* s.Grad(T);
eqs{2} = (1./dt).*((vol-pvMult.*s.pv).*uR-(vol-pvMult0.*s.pv).*uR0) + s.Div( vQ);
% accumulation of energy and advection of enthalpy in fluid, conduction of
% heat is neglected
eqs{2}  =  eqs{2} + ((s.pv/dt).*( pvMult.*uW.*f.rhoWS.*bW - pvMult0.*uW0.*f.rhoWS.*bW0 )...
        +  s.Div( s.faceUpstr(bWvW>0, f.rhoWS.*hW) .* bWvW));
% add energy from pressure boundaries assume outside has sime enthapy as
%
if(~isempty(bcCells))
    qBC=qBC{1};
    hWbc=hW(bcCells);%.*BCtocellMap;
    hWbc(qBC>0)=drivingForces.bc.hW(qBC>0);
    eqs{2} = eqs{2} - BCtocellMap*(hWbc.*f.rhoWS.*qBC);
end
% ensure inflow is correct

if(~isempty(srcCells))
    qSRC=qSRC{1};
    hWsrc=hW(srcCells);
    hWsrc(srcCells(qSRC>0))=drivingForces.src.hW.*f.rhoWs;
    eqs{2}(srcCells) = eqs{2}(srcCells) - qSRC.*hWsrc;
end

% add energy from conduction from boundary
if(~isempty(drivingForces.bcT))
    assert(all(strcmp(drivingForces.bcT.type,'temperature')));
    bc_t=strcmp(drivingForces.bcT.type,'temperature');
    if(any(bc_t))
        %assert(isempty(drivingForces.bc),'can only have temprature boundary with nowflow');
        bc=drivingForces.bcT;
        [bQqQbc,bcT_cell] = temperatureBCContrib(G,s,T,bc);
        eqs{2}(bcT_cell)  = eqs{2}(bcT_cell) + bQqQbc;
    end
end
                                   
                                   
primaryVars = {'pressure', 'T', 'qWs', 'bhp'};
names = {'water', 'temperature'};
types = {'cell', 'cell'};
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
        [cqs, weqs, ctrleqs, wc, state.wellSol, ~]  = wm.computeWellFlux(model, W, wellSol, ...
                                             pBH, {qWs}, pw, rhos, bw, mw, s, {},...
                                             'nonlinearIteration', opt.iteration,'referencePressureIndex', 1);
        
                               
        eqs(3) = weqs;
        eqs{4} = ctrleqs;
        
        eqs{1}(wc) = eqs{1}(wc) - cqs{1};
        % thermal part well
        %%{
        hFw=f.rhoWS*hW(wc);
        hFww=f.rhoWS*wm.Rw*vertcat(W.hW);
        hFw(cqs{1}>0)=hFww(cqs{1}>0);
        %}
        %{
        hFw=f.rhoWS*wm.Rw*vertcat(W.hW);
        hFw(cqs{1}<0)=f.rhoWS*hW(wc(cqs{1}<0));%give wrong derivatives
        %}
        eqs{2}(wc)= eqs{2}(wc) - hFw.*cqs{1};
        names(3:4) = {'waterWells', 'closureWells'};
        types(3:4) = {'perf', 'well'};
    else
        % in reverse mode just gather zero-eqs of correct size
        [eqs(3:4), names(3:4), types(3:4)] = wm.createReverseModeWellEquations(model, state0.wellSol, p0);
        %{
        for eqn = 2:3
            nw = numel(state0.wellSol);
            zw = double2ADI(zeros(nw,1), p0);
            eqs(2:3) = {zw, zw};
        end
        names(2:3) = {'empty', 'empty'};
        types(2:3) = {'none', 'none'};
        %}
    end
end
eqs{2}=eqs{2}/1e6;
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end
%--------------------------------------------------------------------------

function  [bQqQbc,bc_cell] = temperatureBCContrib(G,s,T,bc)        
        assert(all(strcmp(bc.type,'temperature')),'only pressure bc allowed');
        Tr_bc=s.T_r_all(bc.face);
        assert(all(sum(G.faces.neighbors(bc.face,:)>0,2)==1),'bc on internal boundary');
        bc_cell=sum(G.faces.neighbors(bc.face,:),2);        
        Tbc=T(bc_cell);
        dTbc=bc.value-Tbc;
        bQqQbc  = -Tr_bc.*(dTbc);
end









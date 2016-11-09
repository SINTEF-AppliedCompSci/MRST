function [eqs, cq_mass, mix_s, status, cstatus, cq_vol] = computeWellContributionsSingleWell(wellmodel, wellSol, resmodel, q_s, pBH, varw, p, mob, rho, components)

% W = wellmodel.W;
% p = wellmodel.referencePressure;
% b = wellmodel.bfactors;
% r = wellmodel.components;
% m = wellmodel.mobilities;

% perf2well = wellmodel.perf2well;
% Rw = wellmodel.Rw;
% numPh       = numel(b); % # phases
% Tw    = vertcat(W(:).WI);

W = wellmodel.W;
% compi = vertcat(W(:).compi);
assert(numel(wellSol) == 1);
assert(numel(W) == 1);
wc = W.cells;
nc = numel(wc);
numPh = numel(q_s);

% OBS!!
% b = rho;
b = cell(numPh, 1);
rhoS = resmodel.getSurfaceDensities();
for i = 1:numel(b)
    b{i} = rho{i}./rhoS(i);
end

Tw = W.WI;
compi = W.compi;

wellStatus = W.status;
% Perforations/completions are closed if the well are closed or they are
% individually closed
perfStatus = W.cstatus.*wellStatus;
% Closed shut connection by setting WI = 0
Tw(~perfStatus) = 0;

% Well total volume rate at std conds:
qt_s = 0;
for ph = 1:numPh
    qt_s = qt_s + q_s{ph};
end
qt_s = qt_s*wellStatus;

% Get well signs, default should be that wells are not allowed to change sign
% (i.e., prod <-> inj)
if ~wellmodel.allowSignChange % injector <=> w.sign>0, prod <=> w.sign<0
    isInj = W.sign>0;
else
    isInj = double(qt_s)>0;   % sign determined from solution
end

%--------------------------------------------------------------------------
% Pressure drawdown (also used to determine direction of flow)
drawdown    = -(pBH+vertcat(wellSol.cdp)) + p;
connInjInx  = drawdown < 0; %current injecting connections

% A cross-flow connection is is defined as a connection which has opposite
% flow-direction to the well total flow
crossFlowConns = connInjInx ~= isInj;
% If crossflow is not alowed, close connections by setting WI=0
closedConns = ~wellSol.cstatus;
%closedConns    = false(size(crossFlowConns));
if ~wellmodel.allowCrossflow
    closedConns     = or(closedConns, crossFlowConns);
end
Tw(closedConns) = 0;
% Remove closedConns from connInjInx
connInjInx      = and(connInjInx, ~closedConns);

% ------------------ HANDLE FLOW INTO WELLBORE -------------------------
% producing connections phase volumerates:
cq_p = cell(1, numPh);
conEff = ~connInjInx.*Tw;
for ph = 1:numPh
    cq_p{ph} = -conEff.*mob{ph}.*drawdown;
end
% producing connections phase volumerates at standard conditions:
cq_ps = conn2surf(cq_p, b, components, resmodel);
% Sum of phase rates from producing connections at std conds:
q_ps = cell(1, numPh);
for ph = 1:numPh
    q_ps{ph} = sum(cq_ps{ph});
end

isInj = double(qt_s)>0;
% compute avg wellbore phase volumetric rates at std conds.
qt_s_inj = isInj.*qt_s;
wbq = cell(1, numPh);
for ph = 1:numPh
    wbq{ph} = compi(:,ph).*qt_s_inj + ~isInj.*q_s{ph}.*(q_s{ph}>0) - q_ps{ph};
end
% compute wellbore total volumetric rates at std conds.
wbqt = wbq{1};
for ph = 2:numPh
    wbqt = wbqt + wbq{ph};
end
% check for "dead wells":
deadWells = double(wbqt)==0;
if any(deadWells)
    for ph = 1:numPh
        wbq{ph} = wbq{ph}.*(~deadWells) + compi(:, ph).*deadWells;
        % Avoid division by zero
    end
    wbqt(deadWells) = 1;
end
% compute wellbore mixture at std conds
mix_s = cell(1, numPh);
for ph = 1:numPh
    mix_s{ph} = wbq{ph}./wbqt;
end
% ------------------ HANDLE FLOW OUT FROM WELLBORE -----------------------
% total mobilities:
mt = mob{1};
for ph = 2:numPh
    mt = mt + mob{ph};
end
% injecting connections total volumerates
cqt_i = -(connInjInx.*Tw).*(mt.*drawdown);
% volume ratio between connection and standard conditions
volRat  = compVolRat(mix_s, b, components, resmodel);
% injecting connections total volumerates at standard condintions
cqt_is = cqt_i./volRat;
% connection phase volumerates at standard conditions (for output):
cq_s = cell(1,numPh);
for ph = 1:numPh
    cq_s{ph} = cq_ps{ph} + mix_s{ph}.*cqt_is;
end

% Reservoir condition fluxes
cq_vol = cell(1, numPh);
for ph = 1:numPh
    cq_vol{ph} = connInjInx.*cqt_i.*compi(ph) + ~connInjInx.*cq_p{ph};
end
%---------------------- WELL EQUATIONS     -------------------------------
% Well equations
eqs = cell(1, numPh);
for ph = 1:numPh
    eqs{ph} = q_s{ph} - sum(cq_s{ph});
end

if ~all(wellStatus)
    % Overwrite equations with trivial equations for inactive wells
    subs = ~wellStatus;
    for ph = 1:numPh
        eqs{ph}(subs) = q_s{ph}(subs) - double(q_s{ph}(subs));
    end
end
% return mix_s(just values), connection and well status:
mix_s   = cell2mat( cellfun(@double, mix_s, 'UniformOutput', false));
cstatus = ~closedConns;
% For now, don't change status here
status = vertcat(wellSol.status);

cq_mass = cell(1, numPh);
for i = 1:numel(b)
    cq_mass{i} = cq_s{i}.*rhoS(i);
end
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% PRIVATE FUNCTIONS
function qrho = conn2surf(v, b, r, model)
    % in and output cells of ADI
    nPh = numel(v);
    bv = cell(1,nPh);
    for ph = 1:nPh
        bv{ph} = b{ph}.*v{ph};
    end
    qrho = bv;


    dg = isprop(model, 'disgas') && model.disgas;
    vo = isprop(model, 'vapoil') && model.vapoil;
    if (vo || dg)
        [~, isgas] = model.getVariableField('sg');
        [~, isoil] = model.getVariableField('so');
        if isa(model, 'ThreePhaseBlackOilModel')
            if dg
                qrho{isgas} = qrho{isgas} + r{isgas}{isoil}.*bv{isoil};
            end
            if vo
                qrho{isoil} = qrho{isoil} + r{isoil}{isgas}.*bv{isgas};
            end
        elseif dg
            qrho{isgas} = qrho{isgas} + r{isgas}{isoil}.*bv{isoil};
        end
    end
end
%--------------------------------------------------------------------------
function volRat  = compVolRat(cmix_s, b, r, model)
    % first extend mix_s to number of connections:
    nPh = numel(b);
    tmp = cmix_s;


    dg = isprop(model, 'disgas') && model.disgas;
    vo = isprop(model, 'vapoil') && model.vapoil;
    if vo || dg
        [~, isgas] = model.getVariableField('sg');
        [~, isoil] = model.getVariableField('so');
        [rs, rv] = deal(0);
        if dg
            rs = r{isgas}{isoil};
        end
        if vo
            rv = r{isoil}{isgas};
        end
        if (vo || dg) && isa(model, 'ThreePhaseBlackOilModel')
            d = 1-rv.*rs;
            if dg
                tmp{isgas} = (tmp{isgas} - rs.*cmix_s{isoil})./d;
            end
            if vo
                tmp{isoil} = (tmp{isoil} - rv.*cmix_s{isgas})./d;
            end
        elseif dg
            tmp{isgas} = tmp{isgas} - rs.*cmix_s{isoil};
        end
    end

    volRat = 0;
    for ph = 1:nPh
        volRat = volRat + tmp{ph}./b{ph};
    end
end





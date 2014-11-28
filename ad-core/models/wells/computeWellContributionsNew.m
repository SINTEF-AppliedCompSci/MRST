function [eqs, cq_s, mix_s, status, cstatus, Rw, cq_r] = computeWellContributions(wellmodel, model, sol, pBH, q_s)

W = wellmodel.W;
p = wellmodel.referencePressure;
b = wellmodel.bfactors;
r = wellmodel.components;
m = wellmodel.mobilities;

numPh       = numel(b); % # phases
nConn       = cellfun(@numel, {W.cells})'; % # connections of each well
perf2well   = rldecode((1:numel(W))', nConn);
% helpful matrix in dealing with all wells at one go
Rw    = sparse((1:numel(perf2well))', perf2well, 1, numel(perf2well), numel(W));
Tw    = vertcat(W(:).WI);
%active phases
% [~, actPh] = model.getActivePhases();

compi = vertcat(W(:).compi);

cWstatus = vertcat(W(:).cstatus);
% Closed shut connection by setting WI = 0
Tw(~cWstatus) = 0;


% Well total volume rate at std conds:
qt_s = q_s{1};
for ph = 2:numPh
    qt_s = qt_s + q_s{ph};
end
% Get well signs, default should be that wells are not allowed to change sign
% (i.e., prod <-> inj)
if ~wellmodel.allowWellSignChange % injector <=> w.sign>0, prod <=> w.sign<0
    isInj = vertcat(W.sign)>0;
else
    %qt_s =sol.qWs+sol.qOs+sol.qGs;
    isInj = double(qt_s)>0;   % sign determined from solution
end

%--------------------------------------------------------------------------
% Pressure drawdown (also used to determine direction of flow)
drawdown    = -(Rw*pBH+vertcat(sol.cdp)) + p;
connInjInx  = (drawdown <0 ); %current injecting connections

% A cross-flow connection is is defined as a connection which has opposite
% flow-direction to the well total flow
crossFlowConns = connInjInx ~= Rw*isInj;
% If crossflow is not alowed, close connections by setting WI=0
closedConns = ~vertcat(sol.cstatus);
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
for ph = 1:numPh
    cq_p{ph} = -(~connInjInx.*Tw).*(m{ph}.*drawdown);
end
% producing connections phase volumerates at standard conditions:
cq_ps = conn2surf(cq_p, b, r, model);
% Sum of phase rates from producing connections at std conds:
q_ps = cell(1, numPh);
for ph = 1:numPh
    q_ps{ph} = Rw'*cq_ps{ph};
end

isInj = double(qt_s)>0;
% compute avg wellbore phase volumetric rates at std conds.
wbq = cell(1, numPh);
for ph = 1:numPh
    %wbq{ph} = isInj.*q_s{ph} - q_ps{ph};
    wbq{ph} = (isInj.*compi(:,ph)).*qt_s + ~isInj.*q_s{ph}.*(q_s{ph}>0) - q_ps{ph};
%     wbq{ph} = (isInj.*compi(:,ph)).*qt_s - q_ps{ph};
end
% compute wellbore total volumetric rates at std conds.
wbqt = wbq{1};
for ph = 2:numPh
    wbqt = wbqt + wbq{ph};
end
% check for "dead wells":
deadWells =  double(wbqt)==0;
% compute wellbore mixture at std conds
mix_s = cell(1, numPh);
for ph = 1:numPh
    if isa(pBH, 'ADI')
        mix_s{ph} = double2ADI(compi(:,ph), pBH);
    else
        mix_s{ph} = zeros(size(pBH));
    end
    mix_s{ph}(~deadWells) = wbq{ph}(~deadWells)./wbqt(~deadWells);
end
% ------------------ HANDLE FLOW OUT FROM WELLBORE -----------------------
% total mobilities:
mt = m{1};
for ph = 2:numPh
    mt = mt + m{ph};
end
% injecting connections total volumerates
cqt_i = -(connInjInx.*Tw).*(mt.*drawdown);
% volume ratio between connection and standard conditions
volRat  = compVolRat(mix_s, b, r, Rw, model);
% injceting connections total volumerates at standard condintions
cqt_is = cqt_i./volRat;
% connection phase volumerates at standard conditions (for output):
cq_s = cell(1,numPh);
for ph = 1:numPh
    cq_s{ph} = cq_ps{ph} + (Rw*mix_s{ph}).*cqt_is;
end

% Reservoir condition fluxes
cq_r = cell(1, numPh);
for ph = 1:numPh
    cq_r{ph} = connInjInx.*cqt_i.*compi(perf2well,ph) + ~connInjInx.*cq_p{ph};
end
%---------------------- WELL EQUATIONS     -------------------------------
% Well equations
eqs = cell(1, numPh);
for ph = 1:numPh
    eqs{ph} = q_s{ph} - Rw'*cq_s{ph};
end
% return mix_s(just values), connection and well status:
mix_s   = cell2mat( cellfun(@double, mix_s, 'UniformOutput', false));
cstatus = ~closedConns;
%status  = ~deadWells;
% For now, don't change status here
status = vertcat(sol.status);
if(mrstVerbose && any(deadWells) )
%     warning('It exist deadWells')
end
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% PRIVATE FUNCTIONS
function vs = conn2surf(v, b, r, model)
    % in and output cells of ADI
    nPh = numel(v);
    bv = cell(1,nPh);
    for ph = 1:nPh
        bv{ph} = b{ph}.*v{ph};
    end
    vs = bv;


    dg = isprop(model, 'disgas') && model.disgas;
    vo = isprop(model, 'vapoil') && model.vapoil;
    if (vo || dg)
        [~, isgas] = model.getVariableField('sg');
        [~, isoil] = model.getVariableField('so');
        if isa(model, 'ThreePhaseBlackOilModel')
            vs{isgas} = vs{isgas} + r{1}.*bv{isoil};
            vs{isoil} = vs{isoil} + r{2}.*bv{isgas};
        elseif dg
            vs{isgas} = vs{isgas} + r{1}.*bv{isoil};
        end
    end
end
%--------------------------------------------------------------------------
function volRat  = compVolRat(mix_s, b, r, Rw, model)
    % first extend mix_s to number of connections:
    nPh = numel(b);
    cmix_s = cell(1,nPh);
    for ph = 1:nPh
        cmix_s{ph} = Rw*mix_s{ph};
    end
    tmp = cmix_s;


    dg = isprop(model, 'disgas') && model.disgas;
    vo = isprop(model, 'vapoil') && model.vapoil;
    if vo || dg
        [~, isgas] = model.getVariableField('sg');
        [~, isoil] = model.getVariableField('so');
        if (vo || dg) && isa(model, 'ThreePhaseBlackOilModel')
            d = 1-r{1}.*r{2};
            tmp{isgas} = (tmp{isgas} - r{1}.*cmix_s{isoil})./d;
            tmp{isoil} = (tmp{isoil} - r{2}.*cmix_s{isgas})./d;
        elseif dg
            tmp{isgas} = tmp{isgas} - r{1}.*cmix_s{isoil};
        end
    end

    volRat = tmp{1}./b{1};
    for ph = 2:nPh
        volRat = volRat + tmp{ph}./b{ph};
    end
end





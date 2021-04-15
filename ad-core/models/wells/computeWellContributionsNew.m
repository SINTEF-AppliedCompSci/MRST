function [eqs, cq_s, mix_s, status, cstatus, cq_r] = computeWellContributionsNew(wellmodel, model, sol, pBH, q_s)
%Setup well (residual) equations and compute corresponding source terms.
%
% SYNOPSIS:
%   [eqs, cq_s, mix_s, status, cstatus, Rw, cq_r] = ...
%            computeWellContributionsNew(wellmodel, model, sol, pBH, q_s)
%
% PARAMETERS:
%   wellmodel   - Simulation well model.
%   model       - Simulation model.
%   sol         - List of current well solution structures
%   pBH         - Vector of well bhps
%   q_s         - List of vectors of well component volume-rates 
%                 (surface conds) 
%
% RETURNS:
%   eqs         - List of well equations
%   cq_s        - List of vectors containing volumetric component 
%                 source-terms (surface conds). 
%   mix_s       - List of vectors containing volumetric mixture of components 
%                 in wellbroe at connections (surface conds).
%   status      - Logic vector of well statuses
%   cstatus     - Logic vector of well connection statuses
%   Rw          - Sparse matrix representing connection to well mapping
%   cq_r        - List of vectors containing volumetric phase 
%                 source-terms (reservoir conds).
%
% SEE ALSO:
%   WellModel, setupWellControlEquation

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
W = wellmodel.W;
p = wellmodel.referencePressure;
b = wellmodel.bfactors;
r = wellmodel.components;
m = wellmodel.mobilities;

perf2well = wellmodel.perf2well;
Rw = wellmodel.Rw;
numPh       = numel(b); % # phases
Tw    = vertcat(W(:).WI);

compi = vertcat(W(:).compi);

wellStatus = vertcat(W.status);
% Perforations/completions are closed if the well are closed or they are
% individually closed
perfStatus = vertcat(W.cstatus).*wellStatus(perf2well);
% Closed shut connection by setting WI = 0
Tw(~perfStatus) = 0;


% Well total volume rate at std conds:
qt_s = q_s{1}.*wellStatus;
for ph = 2:numPh
    qt_s = qt_s + q_s{ph}.*wellStatus;
end


% Get well signs, default should be that wells are not allowed to change sign
% (i.e., prod <-> inj)
if ~wellmodel.allowWellSignChange % injector <=> w.sign>0, prod <=> w.sign<0
    isInj = vertcat(W.sign)>0;
else
    %qt_s =sol.qWs+sol.qOs+sol.qGs;
    isInj = value(qt_s)>0;   % sign determined from solution
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
conEff = ~connInjInx.*Tw;
for ph = 1:numPh
    cq_p{ph} = -conEff.*m{ph}.*drawdown;
end
% producing connections phase volumerates at standard conditions:
cq_ps = conn2surf(cq_p, b, r, model);
% Sum of phase rates from producing connections at std conds:
q_ps = cell(1, numPh);
for ph = 1:numPh
    q_ps{ph} = Rw'*cq_ps{ph};
end

isInj = value(qt_s)>0;
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
deadWells = value(wbqt)==0;
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
mt = m{1};
for ph = 2:numPh
    mt = mt + m{ph};
end
% injecting connections total volumerates
cqt_i = -(connInjInx.*Tw).*(mt.*drawdown);
% volume ratio between connection and standard conditions
volRat  = compVolRat(mix_s, b, r, Rw, model);
% injecting connections total volumerates at standard condintions
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

if ~all(wellStatus)
    % Overwrite equations with trivial equations for inactive wells
    subs = ~wellStatus;
    for ph = 1:numPh
        eqs{ph}(subs) = q_s{ph}(subs) - value(q_s{ph}(subs));
    end
end
% return mix_s(just values), connection and well status:
mix_s   = cell2mat( cellfun(@value, mix_s, 'UniformOutput', false));
cstatus = ~closedConns;
% For now, don't change status here
status = vertcat(sol.status);
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





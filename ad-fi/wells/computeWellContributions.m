function [eqs, cq_s, mix_s, status, cstatus, Rw] = computeWellContributions(W, sol, pBH, q_s, p, b, r, m, model, ...
                                                             allowWellSignChange, allowCrossflow)
% [eqs, cq_s, closedConns] = getWellContributions(...
%           W, sol, pBH, qt_s, mix_s, p, b, r, m, model, ...
%           allowControlSwitching, allowWellSignChange, allowCrossflow)
% INPUT/OUTPUT
% see getWellContributions

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

numPh       = numel(b); % # phases
nConn       = cellfun(@numel, {W.cells})'; % # connections of each well
perf2well   = rldecode((1:numel(W))', nConn);
% helpful matrix in dealing with all wells at one go
Rw    = sparse((1:numel(perf2well))', perf2well, 1, numel(perf2well), numel(W));
Tw    = vertcat(W(:).WI);
%active phases
actPh = getActivePhases(model);
compi = vertcat(W(:).compi);
compi = compi(:, actPh);

% Get well signs, default should be that wells are not allowed to change sign
% (i.e., prod <-> inj)
if ~allowWellSignChange % injector <=> w.sign>0, prod <=> w.sign<0
    isInj = vertcat(W.sign)>0;
else
   qt_s = q_s{1};
   for ph = 2:numPh
      qt_s = qt_s + q_s{ph};
   end
    isInj = value(qt_s)>0;   % sign determined from solution
end

%--------------------------------------------------------------------------
% Pressure drawdown (also used to determine direction of flow)
drawdown    = -(Rw*pBH+vertcat(sol.cdp)) + p;
connInjInx  = (drawdown <0 ); %current injecting connections

% A cross-flow connection is is defined as a connection which has opposite
% flow-direction to the well total flow
crossFlowConns = connInjInx ~= Rw*isInj;
% If crossflow is not allowed, close connections by setting WI=0
closedConns    = false(size(crossFlowConns));
if ~allowCrossflow
    closedConns     = crossFlowConns;
    Tw(closedConns) = 0;
    % Remove closedConns from connInjInx
    connInjInx      = and(connInjInx, ~closedConns);
end

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
% Well total volume rate at std conds:
qt_s = q_s{1};
for ph = 2:numPh
    qt_s = qt_s + q_s{ph};
end
isInj = value(qt_s)>0;
% compute avg wellbore phase volumetric rates at std conds.
wbq = cell(1, numPh);
for ph = 1:numPh
    %wbq{ph} = isInj.*q_s{ph} - q_ps{ph};
    wbq{ph} = (isInj.*compi(:,ph)).*qt_s - q_ps{ph};
end
% compute wellbore total volumetric rates at std conds.
wbqt = wbq{1};
for ph = 2:numPh
    wbqt = wbqt + wbq{ph};
end
% check for "dead wells":
deadWells =  value(wbqt)==0;
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
% injecting connections total volume rates at standard conditions
cqt_is = cqt_i./volRat;
% connection phase volumerates at standard conditions (for output):
cq_s = cell(1,numPh);
for ph = 1:numPh
    cq_s{ph} = cq_ps{ph} + (Rw*mix_s{ph}).*cqt_is;
end
%---------------------- WELL EQUATIONS     -------------------------------
% Well equations
eqs = cell(1, numPh);
for ph = 1:numPh
    eqs{ph} = q_s{ph} - Rw'*cq_s{ph};
end
% return mix_s(just values), connection and well status:
mix_s   = cell2mat( cellfun(@value, mix_s, 'UniformOutput', false));
cstatus = ~closedConns;
status  = ~deadWells; %any(~deadWells, Rw'*cstatus); % 0 if dead or all conns closed
if mrstVerbose && any(deadWells),
   dW = {W(deadWells).name};
   dW = cellfun(@(w)([w, ' ']), dW, 'UniformOutput', false);
   dW = horzcat(dW{:});
   dW = dW(1:end - 1);
   fprintf('Inactive wells: %s.\n', dW);
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
switch model
    case 'OW' %done
    case 'WG' %done
    case '3P' %done
    case 'BO'
        vs{3} = vs{3} + r{1}.*bv{2};
    case 'VO'
        vs{3} = vs{3} + r{1}.*bv{2};
        vs{2} = vs{2} + r{2}.*bv{3};
    otherwise
        error(['Unknown model: ', model]);
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
switch model
    case 'OW' %done
    case 'WG' %done
    case '3P' %done
    case 'BO'
        tmp{3} = tmp{3} - r{1}.*cmix_s{2};
    case 'VO'
        d = 1-r{1}.*r{2};
        tmp{3} = (tmp{3} - r{1}.*cmix_s{2})./d;
        tmp{2} = (tmp{2} - r{2}.*cmix_s{3})./d;
    otherwise
        error(['Unknown model: ', model]);
end
volRat = tmp{1}./b{1};
for ph = 2:nPh
    volRat = volRat + tmp{ph}./b{ph};
end
end





function [eqs, cq_mass, mix_s, status, cstatus, cq_vol] = computeWellContributionsSingleWell(wellmodel, wellSol, resmodel, q_s, pBH, packed)
% Main internal function for computing well equations and source terms

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

[p, mob, rho, dissolved] = unpackPerforationProperties(packed);

W = wellmodel.W;
assert(numel(wellSol) == 1);
assert(numel(W) == 1);
numPh = numel(q_s);

rhoS = resmodel.getSurfaceDensities();


% Well total volume rate at std conds:
qt_s = 0;
if W.status
    for ph = 1:numPh
        qt_s = qt_s + q_s{ph};
    end
end
drawdown = getDrawdown(wellSol, pBH, p);
[Tw, isInj, connInjInx, cstatus] = getWellTrans(wellmodel, wellSol, qt_s, drawdown);

% ------------------ HANDLE FLOW INTO WELLBORE -------------------------
[cq_vol, cq_s, mix_s_ad] = computePerforationRates(resmodel, wellmodel, connInjInx, q_s, qt_s, drawdown, Tw, isInj, mob, rho, rhoS, dissolved);

%---------------------- WELL EQUATIONS     -------------------------------
% Well equations
eqs = cell(1, numPh);
for ph = 1:numPh
    eqs{ph} = q_s{ph} - sum(cq_s{ph});
end

if ~W.status
    % Overwrite equations with trivial equations for inactive wells
    subs = ~wellStatus;
    for ph = 1:numPh
        eqs{ph}(subs) = q_s{ph}(subs) - value(q_s{ph}(subs));
    end
end
% return mix_s(just values), connection and well status:
mix_s = mix_s_ad;
for i = 1:numel(mix_s)
    mix_s{i} = value(mix_s_ad{i});
end

mix_s   = cell2mat(mix_s);
% For now, don't change status here
status = vertcat(wellSol.status);

cq_mass = cell(1, numPh);
for i = 1:numPh
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
        bv{ph} = reduceToDouble(b{ph}.*v{ph});
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
            d = reduceToDouble(1-rv.*rs);
            if dg
                tmp{isgas} = reduceToDouble(tmp{isgas} - rs.*cmix_s{isoil})./d;
            end
            if vo
                tmp{isoil} = reduceToDouble(tmp{isoil} - rv.*cmix_s{isgas})./d;
            end
        elseif dg
            tmp{isgas} = tmp{isgas} - reduceToDouble(rs.*cmix_s{isoil});
        end
    end

    volRat = 0;
    for ph = 1:nPh
        volRat = volRat + reduceToDouble(tmp{ph}./b{ph});
    end
end

function drawdown = getDrawdown(wellSol, pBH, p)
    % Pressure drawdown (also used to determine direction of flow)
    drawdown    = -(pBH+vertcat(wellSol.cdp)) + p;

end

function [Tw, isInj, connInjInx, cstatus] = getWellTrans(wellmodel, wellSol, qt_s, drawdown)
    W = wellmodel.W;
    Tw = W.WI;

    wellStatus = W.status;
    % Perforations/completions are closed if the well are closed or they are
    % individually closed
    perfStatus = W.cstatus & wellStatus;
    % Closed shut connection by setting WI = 0
    Tw(~perfStatus) = 0;

    % Get well signs, default should be that wells are not allowed to change sign
    % (i.e., prod <-> inj)
    if wellmodel.allowSignChange 
        isInj = value(qt_s)>0;   % sign determined from solution
    else % injector <=> w.sign>0, prod <=> w.sign<0
        isInj = W.sign>0;
    end

    %--------------------------------------------------------------------------
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
    cstatus = ~closedConns;
end

function [cq_vol, cq_s, mix_s] = computePerforationRates(resmodel, wellmodel, connInjInx, q_s, qt_s, drawdown, Tw, isInj, mob, rho, rhoS, dissolved)
    compi = wellmodel.W.compi;

    numPh = numel(rhoS);
    anyInjPerf  = any(connInjInx);
    anyProdPerf = any(~connInjInx);

    b = phaseDensitiesTobfactor(rho, rhoS, dissolved);

    % producing connections phase volumerates:
    if anyProdPerf
        cq_p = cell(1, numPh);

        conEff = ~connInjInx.*Tw;
        for ph = 1:numPh
            cq_p{ph} = -conEff.*mob{ph}.*drawdown;
        end
        % producing connections phase volumerates at standard conditions:
        cq_ps = conn2surf(cq_p, b, dissolved, resmodel);
        
        isInj = value(qt_s)>0;
        wbq = cell(1, numPh);

        for ph = 1:numPh
            wbq{ph} = q_s{ph}.*(q_s{ph}>0) - sum(cq_ps{ph});
        end

        % compute wellbore total volumetric rates at std conds.
        wbqt = 0;
        for ph = 1:numPh
            wbqt = wbqt + wbq{ph};
        end
        % check for "dead wells":
        deadWells = value(wbqt)==0;
        if any(deadWells)
            for ph = 1:numPh
                wbq{ph} = ~deadWells.*wbq{ph} + compi(ph).*deadWells;
                % Avoid division by zero
            end
            wbqt(deadWells) = 1;
        end
        % compute wellbore mixture at std conds
        mix_s = cell(1, numPh);
        for ph = 1:numPh
            mix_s{ph} = wbq{ph}./wbqt;
        end

    else
        mix_s = cell(1, numPh);
        for i = 1:numPh
            mix_s{i} = compi(i);
        end
    end
    
    if anyInjPerf
        % total mobilities:
        mt = 0;
        for ph = 1:numPh
            mt = mt + mob{ph};
        end
        % injecting connections total volumerates
        cqt_i = -(connInjInx.*Tw).*(mt.*drawdown);

        % volume ratio between connection and standard conditions
        volRat  = compVolRat(mix_s, b, dissolved, resmodel);
        % injecting connections total volumerates at standard conditions
        cqt_is = cqt_i./volRat; 
    end
    % connection phase volumerates at standard conditions (for output):
    cq_s = cell(1,numPh);
    % Reservoir condition fluxes
    cq_vol = cell(1, numPh);
    for ph = 1:numPh
        if anyInjPerf && anyProdPerf
            cq_vol{ph} = connInjInx.*cqt_i.*compi(ph) + ~connInjInx.*cq_p{ph};
            cq_s{ph} = cq_ps{ph} + mix_s{ph}.*cqt_is;
        elseif anyInjPerf
            cq_vol{ph} = cqt_i.*compi(ph);
            cq_s{ph} = mix_s{ph}.*cqt_is;
        else
            cq_vol{ph} = cq_p{ph};
            cq_s{ph} = cq_ps{ph};
        end
    end
end


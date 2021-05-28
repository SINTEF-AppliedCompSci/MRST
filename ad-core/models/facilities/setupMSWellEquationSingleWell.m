function [eqs, eqsMS, cq_s, wellSol, mix_s, status, cstatus, cq_r] = setupMSWellEquationSingleWell(wm, model, wellSol0, wellSol, q_s, bhp, pN, alpha, vmix, resProps, dt, iteration)
% Setup well residual equations for multi-segmented wells - and only those.
%
% SYNOPSIS:
%
%   function [eqs, eqsMS, cq_s, mix_s, status, cstatus, cq_r] = setupMSWellEquations(wm, ...
%                                                   model, wellSol0, dt, bhp, q_s, pN, vmix, alpha, wellNo)
%
% PARAMETERS:
%   wm          - Simulation well model.
%   model       - Resservoir Simulation model.
%   wellSol0    - List of well solution structures from previous step
%   dt          - time step
%   bhp         - bottom hole pressure
%   q_s         - volumetric phase injection/production rate
%   pN          - pressure at nodes
%   vmix        - mixture mass flux in segments
%   alpha       - phase composition ratio at nodes
%   wellNo      - Well number identifying the well for which the equations are going to
%                 be assembled (should correspond to a multisegmented well)
%
% RETURNS:
%   eqs         - *Standard* well equations: Cell array of mass balance equations at the connections
%                    - 1 equation for each phase,
%   eqsMS       - multisegment well equations: Cell array of mass balance equations for
%                 each node and of pressure equations
%                    - mass conservation equations for each phase,
%                      Dimension for each equation: number of nodes - 1
%                    - pressure equation
%                      Dimension: number of segments
%   cq_s        - List of vectors containing volumetric component
%                 source-terms (surface conds).
%   mix_s       - List of vectors containing volumetric mixture of components
%                 in wellbores at connections (surface conds).
%   status      - Logic vector of well statuses
%   cstatus     - Logic vector of well connection statuses
%   cq_r        - List of vectors containing volumetric phase
%                 source-terms (reservoir conds).
%
% SEE ALSO:
%   `FacilityModel`

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


    % Get operators
    op = wm.operators;
    w = wm.W;
    
    % properties for connecting cells
    m = resProps.mob;
    b = resProps.b;
    r = resProps.dissolved;
    pr = resProps.pressure;

    % Operators
    numPh = numel(b);
    eqs   = cell(1, numPh);
    eqsMS = cell(1, numPh + 2);

    % Use surface densities as weighting
    rho_s = model.getSurfaceDensities();

    % Compute required properties at nodes
    %[mix_s, rhom] = getNodeProps(model, bhp, q_s, pN, alpha, rho_s);
    [mix_s, mix_r, rhom, mum] = getNodeMix(model, bhp, q_s, pN, alpha, rho_s);
    
    ws0 = wellSol0;
    [bhp0, q_s0, pN0] = deal(ws0.bhp, {ws0.qWs, ws0.qOs, ws0.qGs}, ws0.nodePressure);
    alpha0 = {ws0.nodeComp(:,1), ws0.nodeComp(:,2), ws0.nodeComp(:,3)};
    %[mix_s0, rhom0] = getNodeProps(model, bhp0, q_s0, pN0, alpha0, rho_s);
    [mix_s0, mix_r0, rhom0, mum0] = getNodeMix(model, bhp0, q_s0, pN0, alpha0, rho_s);

    % Pressure drawdown (also used to determine direction of flow)
    corrDP    = 0*(w.cell2node'*rhom).*( (w.connDZ-w.dZ) + (w.cell2node'*(w.nodes.depth - w.refDepth)-w.dZ));
    drawdown  = pr - (w.cell2node'*[bhp; pN])+corrDP ;
    injInx    = (drawdown <0 ); %current injecting connections

    % Connections phase volume rates
    cq = cell(1, numPh);
    Tw = w.WI;

    allowCrossflow = wm.allowCrossflow;
    if ~allowCrossflow
        Tw(injInx) = 0;
    end

    for ph = 1:numPh
        cq{ph} = -Tw.*m{ph}.*drawdown;
    end

    % Connections phase volumerates at standard conditions:
    cq_s = conn2surf(cq, b, r, model);

    % Redefine mobility for injecting connections (totMob*fraction)
    if allowCrossflow
        if any(injInx)
            cq_s_tot = cq_s{1}(injInx);
            for ph = 2:numPh
                cq_s_tot = cq_s_tot + cq_s{ph}(injInx);
            end
            for ph = 1:numPh
                cq_s{ph}(injInx) = cq_s_tot.*(w.cell2node(:, injInx)'*mix_s{ph});
            end
        end
    end

    vols = w.nodes.vol;
    

    up = value(vmix) >= 0;
    % Pressure drop relation
    ddz = op.grad(w.nodes.depth);
    %rhoSeg = rhom(w.segments.topo(:,2));
    % rhoSeg = op.segmentUpstr(up, rhom(2:end));
    rhoSeg = op.aver(rhom);
    muSeg  = op.segmentUpstr(up, mum(2:end));
    dph = norm(gravity())*rhoSeg.*ddz;
    
    for ph = 1:numPh
        ec            = op.div(op.segmentUpstr(up, alpha{ph}).*vmix ) + w.cell2node*cq_s{ph}*rho_s(ph);
        ec(2:end)     = ec(2:end) + (alpha{ph}.*rhom(2:end)-alpha0{ph}.*rhom0(2:end)).*vols(2:end)/dt;
        eqs{ph}       = q_s{ph} - sum(ec)/rho_s(ph);
        eqsMS{ph}     = ec(2:end);
    end
    
    
    % dph = norm(gravity())*op.segmentUpstr(up, rhom(2:end)).*ddz;
    
    % bhp0 = wellSol0(wellNo).bhp;
    if ~isa(w.segments.flowModel, 'function_handle')
        eqsMS{numPh + 1} = (op.grad([bhp; pN]) - dph - wm.pressureDropModel(pN, vmix, alpha, rho_s, rhoSeg))/value(bhp); % Use mixture density
    else % use new form
        eqsMS{numPh + 1} = (op.grad([bhp; pN]) - dph - w.segments.flowModel(vmix, rhoSeg, muSeg))/value(bhp);
    end
    eqsMS{end} = 1;
    for i = 1:numPh
        eqsMS{end} = eqsMS{end} - alpha{i};
    end
    % Update these if necessary
    [status, cstatus, cq_r] = deal(true, true(size(value(pr))), cellfun(@value, cq, 'UniformOutput', false));
    % return mix_s(just values):
    mix_s   = cell2mat( cellfun(@value, mix_s, 'UniformOutput', false));
end

% -------------------------------------------------------------------------

function [b, r, mu, mixr, volRatio] = computeNodeProps(model, p, mixs)
% Assume single PVT-region for now!
    numPh = numel(mixs);
    fluid = model.fluid;
    b  = cell(1,3);
    mu = cell(1,3);
    acp = model.getActivePhases;
    [wat, oil, gas]= deal(acp(1), acp(2), acp(3));

    dg = isprop(model, 'disgas') && model.disgas;
    vo = isprop(model, 'vapoil') && model.vapoil;

    if wat
        [~, wix] = model.getVariableField('sw');
        b{wix}  = fluid.bW(p);
        mu{wix} = fluid.muW(p);
    end

    if oil
        [~, oix] = model.getVariableField('so');
        fo = mixs{oix};
        if ~dg
            b{oix}  = fluid.bO(p);
            mu{oix} = fluid.muO(p);
        end
    end

    if gas
        [~, gix] = model.getVariableField('sg');
        fg = mixs{gix};
        if ~vo
            b{gix}  = fluid.bG(p);
            mu{gix} = fluid.muG(p);
        end
    end

    if dg
        rsMax = fluid.rsSat(p);
        gor = max(abs(fg./fo), 0);
        gor(isnan(value(gor))) = inf;
        rs = min(rsMax, gor);
        b{oix}  = fluid.bO(p, rs, value(rs)>=value(rsMax));
        mu{oix} = fluid.muO(p, rs, value(rs)>=value(rsMax));
    else
        rs = 0;
    end

    if vo
        rvMax = fluid.rvSat(p);
        ogr = max(abs(fo./fg), 0);
        ogr(isnan(value(gor))) = inf;
        rv = min(rvMax, ogr);
        b{gix}  = fluid.bG(p, rv, value(rv)>=value(rvMax));
        mu{gix} = fluid.muG(p, rv, value(rv)>=value(rvMax));
    else
        rv = 0;
    end

    r = {rs, rv};
    
    d = 1-rs.*rv;

    x = mixs;
    if dg
        x{gix} = (x{gix} - rs.*fo)./d;
    end

    if vo
        x{oix} = (x{oix} - rv.*fg)./d;
    end

    volRatio = x{1}./b{1};
    for ph  =2:numPh
        volRatio = volRatio + x{ph}./b{ph};
    end
    if nargout > 1
        mixr = cell(1,3);
        for ph = 1:numPh
            mixr{ph} = x{ph}./(b{ph}.*volRatio);
        end
    end
end

% -------------------------------------------------------------------------

function [mix_s, mix_r, rhom, mum] = getNodeMix(model, bhp, q_s, pN, alpha, rho_s)
    numPh = numel(alpha);
    
    tmp   = alpha;

    % Divide by rho to get volume
    for ph = 1:numel(tmp)
        tmp{ph} = tmp{ph}/rho_s(ph);
    end

    % Total volume
    sumtmp = tmp{1};
    for ph = 2:numPh
        sumtmp = sumtmp + tmp{ph};
    end

    % Volumetric fraction (of component)
    mix_s = cell(1,numPh);
    for ph = 1:numPh
        mix_s{ph} = tmp{ph}./sumtmp;
    end

    % Add top node
    qt_s = abs(q_s{1});
    for ph = 2:numPh, qt_s = qt_s + abs(q_s{ph}); end
    for ph = 1:numPh
        mix_s{ph} = vertcat(abs(q_s{ph})./qt_s, mix_s{ph});
    end

    % compute properties
    [b, r, mu, mix_r, volRatio] = computeNodeProps(model, [bhp; pN], mix_s);
    
    
    % get mixture density
    tmp = mix_s{1}.*rho_s(1);
    for ph = 2:numPh
        tmp = tmp + mix_s{ph}*rho_s(ph);
    end
    rhom     = tmp./volRatio;
    
    % get mixture viscosity
    mum = mix_r{1}.*mu{1};
    for ph = 2:numPh
        mum = mum + mix_r{ph}.*mu{ph};
    end
end


function vs = conn2surf(v, b, r, model)
% In and output cells of ADI
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
            if dg
                vs{isgas} = vs{isgas} + r{isgas}{isoil}.*bv{isoil};
            end
            if vo
                vs{isoil} = vs{isoil} + r{isoil}{isgas}.*bv{isgas};
            end
        elseif dg
            vs{isgas} = vs{isgas} + r{isgas}{isoil}.*bv{isoil};
        end
    end
end

% -------------------------------------------------------------------------


function state = updateStateBlackOilGeneric(model, state, problem, dx, drivingForces)

state0 = state;
W = drivingForces.Wells;
assert(isempty(drivingForces.bc) && isempty(drivingForces.src))



[disgas, vapoil] = deal(false);

if isprop(model, 'vapoil')
    vapoil = model.vapoil;
end

if isprop(model, 'disgas')
    disgas = model.disgas;
end

state = model.updateStateFromIncrement(state, dx, problem,...
                                            'pressure', model.dpMax);

comp = lower(model.componentNames);

satSolVar = intersect(lower(problem.primaryVariables), comp);

if (disgas || vapoil)
    % Black oil with dissolution
    so = model.getProp(state, 'so');
    sw = model.getProp(state, 'sw');
    sg = model.getProp(state, 'sg');
    
    % Magic status flag, see inside for doc
    st = getCellStatus(state0, so, sw, sg, disgas, vapoil);

    dr = model.getIncrement(dx, problem, 'x');
    dsw = model.getIncrement(dx, problem, 'sw');
    % Interpretation of "gas" phase varies from cell to cell, remove
    % everything that isn't sG updates
    dsg = st{3}.*dr - st{2}.*dsw;

    if disgas
        state = model.updateStateFromIncrement(state, st{1}.*dr, problem, 'rs', model.drsMax);
    end
    
    if vapoil
        state = model.updateStateFromIncrement(state, st{2}.*dr, problem, 'rv', model.drsMax);
    end
    
    

%     maxVal = max(abs([dsw, dso, dsg]), [], 2);
%     step   = min(model.dsMax./maxVal, 1);
% 
%     if model.water
%         state = model.setProp(state, 'sw', sw + step.*dsw);
%     end
%     state = model.setProp(state, 'so', so + step.*dso);
%     state = model.setProp(state, 'sg', sg + step.*dsg);
    dso = -(dsg + dsw);
%     if model.water
%         ds = [dsw, dso, dsg];
%     else
%         ds = [dso, dsg];
%     end
    ds = zeros(numel(so), numel(comp));
    ds(:, strcmpi(comp, 'sw')) = dsw;
    ds(:, strcmpi(comp, 'so')) = dso;
    ds(:, strcmpi(comp, 'sg')) = dsg;
    
    state = model.updateStateFromIncrement(state, ds, problem, 's', model.dsMax);
    % We should *NOT* be solving for oil saturation for this to make sense
    assert(~any(strcmpi(satSolVar, 'so')));
    state = computeFlashBlackOil(state, state0, model, st);
    state.s  = bsxfun(@rdivide, state.s, sum(state.s, 2));
else
    % Solution variables should be saturations directly, find the missing
    % link
    fillComponent = setdiff(lower(model.componentNames), satSolVar);
    fillComponent = fillComponent{1};
    % Fill component is whichever component is assumed to fill up the rest of
    % the pores

    solvedFor = ~strcmpi(comp, fillComponent);

    tmp = 1;
    for i = 1:numel(comp)
        if solvedFor(i)
            [state, v] = model.updateStateFromIncrement(state, dx, problem,...
                                                    comp{i}, model.dsMax);
            tmp = tmp - v;
        end
    end
    % Last phase fills the pores. This may make the last phase swell a lot.
    state = model.setProp(state, fillComponent, tmp);
end


% Wells -------------------------------------------------------------------
dqWs = model.getIncrement(dx, problem, 'qWs');
dqOs = model.getIncrement(dx, problem, 'qOs');
dqGs = model.getIncrement(dx, problem, 'qGs');
dpBH = model.getIncrement(dx, problem, 'bhp');

if ~isempty(dpBH)
    dpBH = sign(dpBH).*min(abs(dpBH), abs(model.dpMax.*vertcat(state.wellSol.bhp)));
    
%     wi = strcmpi(comp, 'wg');
%     oi = strcmpi(comp, 'og');
%     gi = strcmpi(comp, 'sg');
    
    for w = 1:numel(state.wellSol)
        ws = state.wellSol(w);
        ws.bhp  = ws.bhp + dpBH(w);
        if model.water
            ws.qWs  = ws.qWs + dqWs(w);
        end
        if model.oil
            ws.qOs  = ws.qOs + dqOs(w);
        end
        if model.gas
            ws.qGs  = ws.qGs + dqGs(w);
        end
        
        tp = ws.type;
        v  = ws.val;
        switch tp
            case 'bhp'
                ws.bhp = v;
            case 'rate'
                % TODO: This uses magic counting and should be fixed, but
                % is dependent on the same being done to the well
                % computations
                ws.qWs = v*W(w).compi(1);
                ws.qOs = v*W(w).compi(2);
                ws.qGs = v*W(w).compi(3);
            case 'orat'
                ws.qOs = v;
            case 'wrat'
                ws.qWs = v;
            case 'grat'
                ws.qGs = v;
        end
        state.wellSol(w) = ws;
    end
end
end

function st = getCellStatus(state, oil, wat, gas, disgas, vapoil)
% Status should be passed on from updateStateVO (to be sure definition is
% identical). rs and rv are assumed to be compatible, i.e. rx = rxSat for
% saturated cells and rx <= rxSat for undersaturated. Three values of
% status are:
% status 0: should not occur (almost water only -> state 3)
% status 1 oil, no gas  : x = rs, sg = 0    , rv = rvMax
% status 2 gas, no oil  : x = rv, sg = 1-sw , rs = rsMax
% status 3 oil and gas  : x = sg, rs = rsMax, rv = rvMax
if isfield(state, 'status')
    status = state.status;
else
    watOnly    = wat > 1- sqrt(eps);
    if ~vapoil
        oilPresent = true;
    else
        oilPresent = or(oil > 0, watOnly);
    end
    if ~disgas
        gasPresent = true;
    else
        gasPresent = or(gas > 0, watOnly);
    end
    status = oilPresent + 2*gasPresent;
end

if ~disgas
    st1 = false;
else
    st1 = status==1;
end
if ~vapoil
    st2 = false;
else
    st2 = status==2;
end
st3 = status == 3;
st = {st1, st2, st3};
end


function wellSol = initWellSolAD(W, model, state0, wellSolInit)
wellSolGiven = (nargin == 4);

if wellSolGiven
    wellSol = wellSolInit;
elseif isfield(state0, 'wellSol')
    wellSol = state0.wellSol;
else
    wellSol = defaultWellSol(state0, W, model);
end
wellSol = assignFromSchedule(W, model, wellSol);
end

function ws = defaultWellSol(state, W, model)
nw = numel(W);
[~, actPh] = model.getActivePhases();

ws = repmat(struct(...
    'name',   [],...
    'status', [],...
    'type',   [],...
    'val',    [],...
    'sign',   [],...
    'bhp',    [],...
    'qWs',    [],...
    'qOs',    [],...
    'qGs',    [],...
    'mixs',   [],...
    'cstatus',[],...
    'cdp',    [],...
    'cqs',    []), [1, nw]);

% just initialize fields that are not assigned in assignFromSchedule
for k = 1:nw
    nConn = numel(W(k).cells);
    nPh   = numel(actPh);
    ws(k).name = W(k).name;
    % To avoid switching off wells, we need to start with a bhp that makes
    % a producer produce and an injector inject. Hence, we intitialize the
    % bhp such that the top connection pressure is 5bar above/below the
    % corresponding well-cell pressure. If W(k).dZ is ~= 0, however, we
    % don't know wht a decent pressure is ...
    % The increment should depend on the problem and the 5bar could be a
    % pit-fall... (also used in initializeBHP in updateConnDP)
    ws(k).bhp = state.pressure(W(k).cells(1)) + 5*W(k).sign*barsa;

    irate = eps;
    if model.water
        ws(k).qWs  = W(k).sign*irate;
    end
    if model.oil
        ws(k).qOs  = W(k).sign*irate;
    end
    if model.gas
        ws(k).qGs  = W(k).sign*irate;
    end
    if isprop(model, 'polymer') && model.polymer
       ws(k).qWPoly = W(k).poly*ws(k).qWs;
    end
    
    ws(k).mixs = W(k).compi(actPh);
    ws(k).qs   = W(k).sign*ones(1, nPh)*irate;
    ws(k).cdp  = zeros(nConn,1);
    ws(k).cqs  = zeros(nConn,nPh);
    % Additional model dependent fields
    if isprop(model, 'polymer') && model.polymer% polymer model
       ws(k).qWPoly = 0;
    end
end
end

function ws = assignFromSchedule(W, model, ws)
% set fields that should be updated if control has changed
for k = 1:numel(W)
    ws(k).status  = W(k).status;
    ws(k).type    = W(k).type;
    ws(k).val     = W(k).val;
    ws(k).sign    = W(k).sign;
    ws(k).cstatus = W(k).cstatus;

    tp = W(k).type;
    v  = W(k).val;

    switch tp
        case 'bhp'
            ws(k).bhp = v;
        case 'rate'
            if model.water
                ws(k).qWs = v*W(k).compi(1);
            end
            if model.oil
                ws(k).qOs = v*W(k).compi(2);
            end
            if model.gas
                ws(k).qGs = v*W(k).compi(3);
            end
            if isprop(model, 'polymer') && model.polymer
                ws(k).qWPoly = W(k).poly*ws(k).qWs;
            end
        case 'orat'
            ws(k).qOs = v;
        case 'wrat'
            ws(k).qWs = v;
        case 'grat'
            ws(k).qGs = v;
    end % No good guess for qOs, etc...
end
end


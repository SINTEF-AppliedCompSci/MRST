function [state, varargout] = solveStationaryPressureStein(G, state, system, W, mob, pv, T, N, varargin)

opt = struct('objective',       [],...
             'maxiter',         10,...
             'src',             zeros(G.cells.num, 1),...
             'linsolve',        @mldivide,...
             'computeTracer',   true, ...
             'subset',          []);

opt = merge_options(opt, varargin{:});
if isempty(opt.subset)
    inx = ':';
else
    inx = opt.subset;
end

s = system.s;

converged = false;
j = 0;
perf2well = wellPerf(W);

%if ~isfield(state, 'wellSol')
    state.wellSol = initWellSol(W, mean(state.pressure));
%end

%while ~converged && j < opt.maxiter
    %disp('*')
    bhp = vertcat(state.wellSol.pressure);
    fluxes = vertcat(state.wellSol.flux);


    [pressure, fluxes, bhp]  = initVariablesADI(state.pressure, fluxes, bhp);

    % Solve pressure
    [eq_p faceMob] = assemblePressureEq(state, G, W, T, N, pressure, mob, s, fluxes, bhp, pv, opt, [], []);


    converged = all(cellfun(@(x) norm(x.val, inf), eq_p) < 1e-5);
%     converged = norm(eq_p{1}.val, inf)/norm(eq_p{1}.jac{1}, inf) < 1e-8;
%     if converged && j > 1; break; end;
%    if j > 1; break; end

    sol = solveEqs(eq_p, opt.linsolve);

    j = j + 1;
    % Update state
    pressure = pressure + sol{1};
    state.pressure = pressure.val;
    state.flux = faceMob.*s.grad(pressure.val);

    %state.flux(~any(G.faces.neighbors == 0, 2)) = flux.val;
    %state.flux = reshape(state.flux, [], 1);
    for i = 1:numel(W)
        state.wellSol(i).flux = state.wellSol(i).flux + sol{2}(perf2well == i);
        state.wellSol(i).pressure = state.wellSol(i).pressure + sol{3}(i);
    end
%end
    bhp = vertcat(state.wellSol.pressure);
    fluxes = vertcat(state.wellSol.flux);
if j > 1
%     warning('mrst:incompPressureIter', 'Convergence should not require iterations for the present model! (Used %d)', j)
end



if nargout > 1
    % SOLVE TOF
    % make inx
    isPos = (state.flux>0);
    NN = N;
    NN(~isPos, :) = N(~isPos, [2,1]);
    th = max(abs(state.flux))*1e-7;
    M  = sparse(NN(:,2), NN(:,1), double(abs(state.flux)>th), G.cells.num, G.cells.num);
    inx = false(G.cells.num,1);
    inx(vertcat(W.cells)) = true;
    ninx = 0;
    while ninx ~= nnz(inx)
        ninx = nnz(inx);
        inx = logical(inx + M*inx);
    end


%    [pressure, fluxes, bhp, tau_forward, tau_backward]  = initVariablesADI(pressure.val, fluxes.val, bhp.val, ones(G.cells.num, 1), ones(G.cells.num, 1));
    [pressure, fluxes, bhp, tau_forward, tau_backward]  = initVariablesADI(state.pressure, fluxes, bhp, zeros(G.cells.num, 1), zeros(G.cells.num, 1));
    [eq, faceMob] = assemblePressureEq(state, G, W, T, N, pressure, mob, s, fluxes, bhp, pv, opt, tau_forward, tau_backward);


    D = SolveTOFEqsADIStein(eq, state, W, opt.computeTracer, inx);
    varargout{1} = D;
end

if nargout > 2


    assert(~isempty(opt.objective));
    scaling.well = getWellScaling(W);
    % equations must be updated since tof has been updated
    [pressure, fluxes, bhp, tau_forward, tau_backward]  = initVariablesADI(state.pressure, double(fluxes), double(bhp), D.tof(:,1), D.tof(:,2));
    [eq, faceMob] = assemblePressureEq(state, G, W, T, N, pressure, mob, s, fluxes, bhp, pv, opt, tau_forward, tau_backward);
    varargout{2} = SolveAdjointTOFEqs(eq, D, opt.objective(state, D), scaling, perf2well, inx);
end

end



function [eq faceTransMob] = assemblePressureEq(state, G, W, T, N, pressure, totMob, s, fluxes, bhp, pv, opt, tau_forward, tau_backward)

perf2well = wellPerf(W);

findTof = ~isempty(tau_forward);


%totMob = getTotalMobility(mob, state, pressure);
%innerf = ~any(G.faces.neighbors == 0, 2);

%fn = G.faces.neighbors(innerf, :);


%f2hf = face2halfface(G);
%f2hf = f2hf(innerf, :);

%tm1 = totMob(fn(:,1)).*T(f2hf(:, 1));
%tm2 = totMob(fn(:,2)).*T(f2hf(:, 2));

faceTransMob = (2 ./ (1./totMob(N(:,1)) + 1./totMob(N(:,2)))).*T;


flux = faceTransMob.*s.grad(pressure);
pressureeq = s.div(flux) - opt.src;

if findTof
    forward_tof  = s.div(flux.*s.faceUpstr(s.grad(pressure) >= 0, tau_forward)) - pv;
    backward_tof = s.div(-flux.*s.faceUpstr(s.grad(pressure) < 0, tau_backward)) - pv;
    isInj = arrayfun(@(x) sum(x.flux) >= 0, state.wellSol) .';
    isProd = ~isInj;
end


well_closure = 0*bhp;
peaceman = 0*fluxes;
wscale = getWellScaling(W);

for i = 1:numel(W)
    w = W(i);
    wc = w.cells;
    wperf = find(perf2well == i);
    wflux = fluxes(wperf);

    if strcmpi(w.type, 'bhp')
        well_closure(i) = w.val - bhp(i);
    else
        % Ensure that sum of perforation fluxes is equal to prescribed rate
        % for each well
        well_closure(i) = w.val;
        for j = 1:numel(wperf)
            well_closure(i) = well_closure(i) - wflux(j);
        end
    end
    %well_closure(i) = well_closure(i)*wscale(i);

    % Subtract source terms from pressure equation
    pressureeq(w.cells) = pressureeq(w.cells) - wflux;

    % Assemble peaceman model for each perforation
    peaceman(wperf) = wflux - w.WI.*totMob(wc).*(bhp(i) - pressure(wc));

    if findTof
        if isProd(i)
            forward_tof(wc)  = forward_tof(wc) - wflux.*tau_forward(wc);
        elseif isInj(i)
            backward_tof(wc) = backward_tof(wc) + wflux.*tau_backward(wc);
        end
    end
end

if ~any(strcmpi({W.type}, 'bhp'))
   % if mrstVerbose
   %      warning('mrst:badpressure', 'No bottom hole pressure controls - pressure is not well defined');
   % end
   pressureeq(1) = pressureeq(1) + pressure(1);
end

if findTof
    eq = {pressureeq, peaceman, well_closure, forward_tof, backward_tof};
else
    eq = {pressureeq, peaceman, well_closure};
end

end

function totMob = getTotalMobility(fluid, state, pressure)
    sW   = state.s(:,1);
    if size(state.s, 2) > 2
        sG   = state.s(:,3);
    else
        sG = zeros(size(state, 1), 1);
    end
    if isfield(state, 'rs')
        rs   = state.rs;
    else
        rs = zeros(size(state.pressure));
    end
    isSat = (sG>0) | (1 - sW + sG)  == 0;

    if abs(nargin(fluid.relPerm)) == 2
        [krW, krO] = fluid.relPerm(sW);
        krG = 0;
    else
        [krW, krO, krG] = fluid.relPerm(sW, sG);
    end

    l_wat = krW./fluid.muW(pressure);
    l_oil = krO./fluid.muO(pressure, rs, isSat);
    l_gas = krG./fluid.muG(pressure);
    totMob = l_oil + l_wat + l_gas;
end

function dx = solveEqs(eqs, linsolve)
    %eqs{1} = eqs{1}(inx1);
    %eqs{1}.jac{1} = eqs{1}.jac{1}(:, inx1);

    numVars = cellfun(@numval, eqs)';
    cumVars = cumsum(numVars);
    ii = [[1;cumVars(1:end-1)+1], cumVars];

    %eliminate rates
    [eqs, eq_r] = elimVars(eqs, 2);
    %eliminate closure
    [eqs, eq_w] = elimVars(eqs, 2);

    eqs_c = cat(eqs{:});

    J = -eqs_c.jac{:};
    %tmp      = zeros(ii(end,end), 1);
    tmp = linsolve(J, eqs_c.val);

    dx{1} = tmp;

    %recover
    dx{3} = recoverVars(eq_w, 2,    {dx{1}});
    dx{2} = recoverVars(eq_r, 2,    {dx{1}, dx{3}});

    %eqn = size(ii,1);
    %dx = cell(eqn,1);
    %for i = 1:eqn
    %        dx{i} = tmp(ii(i,1):ii(i,2));
    %end
end

function grad = SolveAdjointTOFEqs(eqs, D, objk, scaling, perf2well, inx)
    ni = size(D.itracer, 2);
    np = size(D.ptracer, 2);
    % reduce eqs acc to inx:
    eqs{1} = eqs{1}(inx);
    eqs{4} = eqs{4}(inx);
    eqs{5} = eqs{5}(inx);
    for k = 1:numel(eqs)
        eqs{k}.jac{1} = eqs{k}.jac{1}(:,inx);
        eqs{k}.jac{4} = eqs{k}.jac{4}(:,inx);
        eqs{k}.jac{5} = eqs{k}.jac{5}(:,inx);
    end

    for k = 1:numel(objk.jac)
        if ~ismember(k, [2,3]);
            objk.jac{k} = objk.jac{k}(:,inx);
        end
    end
    numVars = cellfun(@numval, eqs)';

    % First three equations are pressure, well flux and well closure
    % Fourth is forward TOF
    % Fifth is backward tof
    % 4 + ni is injection tracers
    % 5 + ni + 1 is backward tof
    % 5 + ni + 1 -> end is production tracers
    i_f = 4;
    %l_forward = tofRobustFix(-eqs{4}.jac{4}) .' \ flattenVector(objk.jac, i_f:i_f+ni);
    l_forward = -eqs{4}.jac{4} .' \ flattenVector(objk.jac, i_f:i_f+ni);

    l_forward(~isfinite(l_forward)) = deal(0);

    i_b = 4 + ni + 1;
    %l_backward = tofRobustFix(-eqs{5}.jac{5}) .' \flattenVector(objk.jac, i_b:i_b+np);
    l_backward = -eqs{5}.jac{5} .' \flattenVector(objk.jac, i_b:i_b+np);

    l_backward(~isfinite(l_backward)) = deal(0);

    p_ind = 1:3;
    A = flattenJacobian(eqs{1}.jac, p_ind);
    q = flattenJacobian(eqs{2}.jac, p_ind);
    c = flattenJacobian(eqs{3}.jac, p_ind);

    A_sys = [A; q; c] .';

    rhs_sys = full(flattenJacobian(objk.jac, p_ind)) .';

    q_ind = numVars(1) + 1 : sum(numVars(1:2));
    pressure_i = 1:numVars(1);

    % Add in dTOF/dPerfFlux terms to right hand side
    rhs_sys(q_ind) = rhs_sys(q_ind) + sum(eqs{4}.jac{2}.' * l_forward, 2)  + sum(eqs{5}.jac{2}.' * l_backward, 2);
    % dTOF/dPressure terms...
    rhs_sys(pressure_i) = rhs_sys(pressure_i) + sum(eqs{4}.jac{1} .' * l_forward, 2)...
                                              + sum(eqs{5}.jac{1} .' * l_backward, 2) ;


    l_p = -A_sys\rhs_sys;

    grad.pressure = l_p(1:numVars(1));
    grad.fluxes = l_p(q_ind);
    %grad.well = scaling.well.*l_p(sum(numVars(1:2)) + 1: sum(numVars(1:3)));
    grad.well = l_p(sum(numVars(1:2)) + 1: sum(numVars(1:3)));
    grad.tof = [l_forward(:, 1), l_backward(:, 1)];
    grad.itracer = l_forward(:, 2:end);
    grad.ptracer = l_backward(:, 2:end);
    grad.objective = objk;
end

function J = flattenJacobian(J, index)
    J = J(index);
    J = horzcat(J{:});
end

function v = flattenVector(v, index)
    v = v(index);
    v = full(vertcat(v{:})) .';
end

function s = getWellScaling(W)
    s = ones(numel(W), 1);
    s(strcmpi({W.type}, 'bhp')) = mean(vertcat(W.WI));
end

function perf2well = wellPerf(W)
    nPerf  = arrayfun(@(x)numel(x.cells), W)';
    perf2well = rldecode((1:numel(W))', nPerf);
end


function f2hf = face2halfface(G)
nf     = diff(G.cells.facePos);
cellno = rldecode(1 : G.cells.num, nf, 2) .';
t      = G.faces.neighbors(G.cells.faces(:,1), 1) == cellno;
f2hf   = accumarray([double(G.cells.faces(:,1)), double(2 - t)], ...
    (1 : numel(cellno)) .', [G.faces.num, 2]);

end

%--------------------------------------------------------------------------
function [eqs, eqn] = elimVars(eqs, n)
% eliminate set of unknowns nr n using equation n ()
solveInx = setdiff(1:numel(eqs), n);
eqn      = eqs{n};

for eqNum = solveInx
    for jacNum = solveInx
        eqs{eqNum}.jac{jacNum} = eqs{eqNum}.jac{jacNum} - eqs{eqNum}.jac{n}*(eqn.jac{n}\eqn.jac{jacNum});
    end
    eqs{eqNum}.val = eqs{eqNum}.val - eqs{eqNum}.jac{n}*(eqn.jac{n}\eqn.val);
end

eqs  = eqs(solveInx);
for eqNum = 1:numel(eqs)
    eqs{eqNum}.jac = eqs{eqNum}.jac(solveInx);
end

end
%--------------------------------------------------------------------------
function x = recoverVars(eq, n, sol)
% recover variables x at position n using solutions sol
solInx = [1:(n-1), (n+1):(numel(sol)+1)];
x = - eq.jac{n}\(eq.val);
for k  = 1:numel(solInx)
    x = x - eq.jac{n}\(eq.jac{solInx(k)}*sol{k});
end
end

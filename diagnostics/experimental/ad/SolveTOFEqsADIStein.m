function D = SolveTOFEqsADIStein(eqs, state, W, computeTracer, varargin)

[tof1 itracer inj] = solveTOF(eqs, 'forward', state, W, computeTracer, varargin{:});
[tof2 ptracer prod] = solveTOF(eqs, 'backward', state, W, computeTracer, varargin{:});

D.tof = [tof1, tof2];
D.itracer = itracer;
D.ptracer = ptracer;
D.inj = inj;
D.prod = prod;

if computeTracer
    [~,I]      = sort(D.ptracer,2,'descend');
    D.ppart    = I(:,1);

    [~,I]     = sort(D.itracer,2,'descend');
    D.ipart   = I(:,1);
else
    D.ppart = [];
    D.ipart = [];
end

end



function [tof, tracer, isInj] = solveTOF(eqs, direction, state, W, computeTracer, varargin)

    if strcmpi(direction, 'forward')
        self = 4;
        other = 5;
        sgn = 1;
    else
        self = 5;
        other = 4;
        sgn = -1;
    end

    if nargin > 5
        inx = varargin{:};
    else
        inx = (1:numel(state.pressure))';
    end
    isInj = find(arrayfun(@(x) sgn*sum(x.flux) >= 0, state.wellSol)) .';
    inxMap = zeros(numel(state.pressure),1);
    inxMap(inx) = (1:nnz(inx))';
    e = eqs{self}(inx);
    rhs = e.val;
    %A = tofRobustFix(-e.jac{self});
    A = -e.jac{self}(:,inx);
    %A = A(inx, inx);
    % Well flux and pressure influence on local system
    %adjust = e.jac{2}*sgn*vertcat(state.wellSol.flux) - e.jac{1}(:,inx)*state.pressure(inx);
    %adjust = -e.jac{2}*vertcat(state.wellSol.flux) - e.jac{1}(:,inx)*state.pressure(inx);

    % These are zero as long as initial guess is zero !!!
    adjust = e.jac{2}*vertcat(state.wellSol.flux) + e.jac{1}(:,inx)*state.pressure(inx);

    % Well flux influence on the other time of flight system. We will use
    % the inverted source terms from the other system to drive tracer flow.
    adjust_tracer = eqs{other}.jac{2}(inx,:)*sgn*vertcat(state.wellSol.flux);

    rhs = rhs + adjust;
    if computeTracer
        rhs = [rhs, zeros(numel(rhs), numel(isInj))];
        for i = 1:numel(isInj)
            wc = inxMap(W(isInj(i)).cells);
            tmp = 0*e.val;
            tmp(wc) = -sgn;

            tmp = tmp.*adjust_tracer;

            rhs(:, i+1) = tmp;
        end
    end
    sol = zeros(numel(state.pressure), size(rhs,2));
    sol(:,1) = 1000*year;
    sol(inx,:) = A\rhs;
    %sol(isnan(sol)) = deal(0);
    %sol(isinf(sol)) = deal(max(sol(~isinf(sol))));

    tof = sol(:,1);
    tof(tof > 1000*year) = 1000*year;
    if computeTracer
        %tracer = min(max(sol(:, 2:end), 0), 1);
        tracer = sol(:, 2:end);
    else
        tracer = [];
    end
end

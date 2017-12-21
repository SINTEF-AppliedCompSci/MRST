classdef PressureReducedLinearSystem < ReducedLinearizedSystem
    properties
        accumulationTerms
        flashEquations
        equations0
        assembled
        wellVarIndices
        nwellvar
        wellvars
        wellvarNames
        w
        model
    end
    
    methods
        function problem = PressureReducedLinearSystem(varargin)
            problem = problem@ReducedLinearizedSystem(varargin{:});
            problem.assembled = false;
        end
        
        function problem = assembleSystem(problem)
            if ~problem.assembled
                problem = assembleSystem@ReducedLinearizedSystem(problem);
                problem.w = problem.getWeights();
%                 problem.w = getImpesWeightsOverallComposition(problem.model, problem.state, 1);
                
                weights = problem.w;
                
                Ap = sparse(0);
                
                A0 = problem.A;
                b0 = problem.b;
                
                q = 0;
                [ncell, ncomp] = size(weights);
                ex = [1:ncell, problem.wellVarIndices];
%                 ex = 1:ncell;

                eqs = cell(ncomp, 1);
                W = cell(ncomp, 1);
                
                
                eq = 0;
                nder = numel(ex);
                if ~isempty(A0)
                    for i = 1:ncomp
                        ix = (1:ncell) + (i-1)*ncell;
                        eqs{i} = ADI(b0(ix), -A0(ix, ex));

                        if problem.iterationNo == 1 || ~isfield(problem.state, 'w_p')
                            wder = sparse([], [], [], ncell, nder);
                        else
                            dp = problem.state.pressure - problem.state.w_p;
                            ddp = (weights(:, i) - problem.state.w(:, i))./dp;
                            ddp(~isfinite(ddp)) = 0;
                            wder = sparse(1:ncell, 1:ncell, ddp, ncell, nder);
                        end

                        W{i} = ADI(weights(:, i), wder);
                        eq = eq + W{i}.*eqs{i};
                    end
                    problem.b = eq.val;
                    Ap = -eq.jac{1};
                else
                    problem.b = 0;
                    for i = 1:ncomp
                        ix = (1:ncell) + (i-1)*ncell;
                        problem.b = problem.b + weights(:, i).*b0(ix);
                    end
                end
                problem.b = [problem.b; b0(problem.wellVarIndices)];
                problem.assembled = true;
                if isempty(problem.A)
                    % res only mode
                    return
                end
                problem.A = [Ap; A0(problem.wellVarIndices, ex)];
                problem.D = [];
                problem.E = [];
                problem.h = [];
            end
        end

        function [dp, report] = processResultAfterSolve(problem, dp, report)
            [dp, report] = processResultAfterSolve@ReducedLinearizedSystem(problem, dp, report);
            m = problem.model;
            nc = m.G.cells.num;
            if 1
                dp(1:nc) = m.limitUpdateRelative(dp(1:nc), problem.state.pressure, m.dpMaxRel);
            else
                [~, rcp] = m.limitUpdateRelative(dp(1:nc), problem.state.pressure, m.dpMaxRel);
                dp = dp*min(rcp);
            end
            
            [dx, dy, ds, dL, twoPhase] = getUpdatesFromPressure(m, problem.state, dp);
%             dp(twoPhase) = dp(twoPhase).*wcomp;
            
            w_var = cellfun(@(x) ['w_', x], m.EOSModel.fluid.names, 'uniformoutput', false);
            v_var = cellfun(@(x) ['v_', x], m.EOSModel.fluid.names, 'uniformoutput', false);
            assert(strcmp(problem.primaryVariables{1}, 'pressure'))
            
            
            ix = find(cellfun(@(x) isa(x, 'ADI'), problem.equations), 1);
            % Calculate positions in newton increment
            numVars = problem.equations{ix}.getNumVars();
            dp0 = dp;
           
            nc = m.G.cells.num;
            well_dp = dp(nc+1:end);
            dpressure = dp(1:nc);
            
            cvars = cumsum([1; numVars]);
            dp = zeros(sum(numVars), 1);
            
            for i = 1:numel(problem.primaryVariables)
                var = problem.primaryVariables{i};
                ix = cvars(i):cvars(i+1)-1;
                switch var
                    case 'pressure'
                        next = dpressure;
                    case 'sato'
                        next = ds(:, 1 + m.water);
                    case 'satg'
                        next = ds(:, 2 + m.water);
                    case 'satw'
                        next = zeros(nc, 1);
                        next(twoPhase) = ds(:, 1);
                    case w_var
                        if any(twoPhase)
                            next = dy(:, strcmp(var, w_var));
                        else
                            next = [];
                        end
                    case v_var
                        next = zeros(nc, 1);
                        aa = strcmp(var, v_var);
                        next(twoPhase) = dx(:, aa);
                    otherwise
                        continue
                end
                dp(ix) = next;
            end
            if 0
                dp_up = dp;
                dp_up(1:numel(dp0)) = 0;
    %             dp_up((nc+1):numel(dp0)) = 0;
                offset = 0;
                for i = 1:numel(problem.wellvarNames)
                    act = find(strcmpi(problem.primaryVariables, problem.wellvarNames{i}));
                    subs = cvars(act):cvars(act+1)-1;
                    v = well_dp(offset + (1:numel(subs)));

                    offset = offset + numel(subs);

                    var = problem.wellvars{i};
                    var = combineEquations(var);

                    v = v + var.jac{1}*dp_up;

                    dp(subs) = v;
                end
            else
                % Re-compute well equations, accounting for all changes in
                % saturation, composition, etc
                isWell = strcmpi(problem.types, 'perf') | strcmpi(problem.types, 'well');
                if any(isWell)
                    wEqs = combineEquations(problem.equations{isWell});
                    Aw = -wEqs.jac{1};
                    dp_up = dp;
                    dp_up((nc+1):numel(dp0)) = 0;

                    q = wEqs.val - Aw*dp_up;
                    Aw = Aw(:, problem.wellVarIndices);

                    d_well = Aw\q;
                    dp(problem.wellVarIndices) = d_well;
                end
            end
        end

        function w = getWeights(problem)
            acc = problem.accumulationTerms;
            state = problem.state;
            [ncell, ncomp] = size(problem.state.components);
            hasWater = size(state.s, 2) == 3;
            if hasWater
                ncomp = ncomp + 1;
            end
            c = combineEquations(acc{:});
            if isnumeric(c)
                w = ones(ncell, ncomp);
                return;
            end
            J = c.jac{1};

            if ~isempty(problem.reorder)
                J = J(:, problem.reorder);
            end
            
            J(:, problem.wellVarIndices) = [];
            
            ndof = ncell*ncomp;
            [B, C, D, E] = getBlocks(J, ndof);
            [L, U] = lu(E);
            A = B - C*(U\(L\D));
%             A = B - C*(E\D);
            b = zeros(ndof, 1);
            b(1:ncell) = 1/barsa;

            w = (A')\b;
            w = reshape(w, [], ncomp);
            w = bsxfun(@rdivide, w, sum(abs(w), 2));
            w = bsxfun(@rdivide, w, sum(state.rho.*state.s, 2));
        end
    end
end

function [B, C, D, E] = getBlocks(J, ndof)
    start = 1:ndof;

    if 0
        stop = (ndof+1):size(J, 2);
        B = J(start, start);
        C = J(start, stop);
        D = J(stop, start);
        E = J(stop, stop);
    else
       [ix, jx, vx] = find(J);
       n = size(J, 2);
       keep = false(n, 1);
       keep(start) = true;
       nk = ndof;

       keepRow = keep(ix);
       keepCol = keep(jx);
       kb = keepRow & keepCol;
       B = sparse(ix(kb), jx(kb), vx(kb), nk, nk);

       kc = keepRow & ~keepCol;
       C = sparse(ix(kc), jx(kc) - nk, vx(kc), nk, n - nk);

       kd = ~keepRow & keepCol;
       D = sparse(ix(kd) - nk, jx(kd), vx(kd), n - nk, nk);

       ke = ~keepRow & ~keepCol;
       E = sparse(ix(ke) - nk, jx(ke) - nk, vx(ke), n - nk, n - nk);
    end

end

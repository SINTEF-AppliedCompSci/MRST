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
                pmodel = problem.model;
                state = problem.state;
                acc = problem.accumulationTerms;
                [weights, dwdp] = getPartialVolumes(pmodel, state, acc, ...
                    'reorder',                    problem.reorder, ...
                    'wellVarIndices',             problem.wellVarIndices, ...
                    'iteration',                  problem.iterationNo, ...
                    'singlePhaseStrategy',        pmodel.singlePhaseStrategy, ...
                    'twoPhaseStrategy',           pmodel.twoPhaseStrategy, ...
                    'singlePhaseDifferentiation', pmodel.singlePhaseDifferentiation, ...
                    'twoPhaseDifferentiation',    pmodel.twoPhaseDifferentiation);
                
                problem.w = weights;
                Ap = sparse(0);
                
                A0 = problem.A;
                b0 = problem.b;
                
                q = 0;
                [ncell, ncomp] = size(weights);
                ex = [1:ncell, problem.wellVarIndices];

                eqs = cell(ncomp, 1);
                W = cell(ncomp, 1);
                eq = 0;
                nder = numel(ex);
                if ~isempty(A0)
                    for i = 1:ncomp
                        ix = (1:ncell) + (i-1)*ncell;
                        eqs{i} = ADI(b0(ix), -A0(ix, ex));

                        if isempty(dwdp)
                            wder = sparse([], [], [], ncell, nder);
                        else
                            wder = sparse(1:ncell, 1:ncell, dwdp(:, i), ncell, nder);
                        end
                        wi = weights(:, i);

                        W{i} = ADI(wi, wder);
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
            
            cnames = m.EOSModel.getComponentNames();
            w_var = cellfun(@(x) ['w_', x], cnames, 'uniformoutput', false);
            v_var = cellfun(@(x) ['v_', x], cnames, 'uniformoutput', false);
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

                    q = wEqs.val - Aw*dp_up;
                    Aw = Aw(:, problem.wellVarIndices);

                    d_well = Aw\q;
                    dp(problem.wellVarIndices) = d_well;
                end
            end
        end

    end
end

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

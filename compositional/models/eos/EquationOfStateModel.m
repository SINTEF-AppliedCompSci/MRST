classdef EquationOfStateModel < PhysicalModel
    % Equation of state model. Implements generalized two-parameter cubic
    % equation of state with Newton and successive substitution solvers, as
    % well as standard functions for computing density and viscosity.
    properties
        CompositionalMixture % CompositionalMixture
        omegaA % Parameter for EOS
        omegaB % Parameter for EOS
        method = 'ssi' % Type of method
        maxSSI = inf;
        PropertyModel % Model to be used for property evaluations
        selectGibbsMinimum = true; % Use minimum Gibbs energy to select Z
        alpha = [];
        minimumComposition = 1e-8; % Minimum composition value (for numerical stability)
        minimumSaturation  = 1e-8; % Minimum total EOS phase saturation 
        equilibriumConstantFunctions = {};
        extraOutput = 0;
    end

    properties (Access = private)
        eosA
        eosB
        eosType
    end

    methods
        function model = EquationOfStateModel(G, fluid, eosname)
            if nargin < 3
                eosname = 'PR';
            end
            model = model@PhysicalModel(G);
            assert(isa(fluid, 'CompositionalMixture'));
            model.CompositionalMixture = fluid;
            model.nonlinearTolerance = 1e-4;
            model.PropertyModel = CompositionalPropertyModel(fluid);
            model = model.setType(eosname);
        end
        
        function n = getNumberOfComponents(model)
            n = model.CompositionalMixture.getNumberOfComponents();
        end

        function n = getComponentNames(model)
            n = model.CompositionalMixture.names;
        end
        
        function Z = solveCubicEOS(model, A, B)
            % Peng Robinson equation of state in form used by Coats
            [E2, E1, E0] = model.getCubicCoefficients(A, B);
            if 1
                % Use vectorized cubic solver
                Z = cubicPositive(E2, E1, E0);
            else
                % Fallback to Matlab's root solver
                n = numel(A);
                Z = zeros(n, 3);
                for i = 1:n
                    Z(i, :) = roots([1, E2(i), E1(i), E0(i)]);
                end
            end
            % Avoid picking any complex roots. At least one root will be
            % real. We replace these invalid values with NaN so that when
            % min or max is used, they will ignore any previously complex
            % values.
            mtol = 0;
            Z(abs(imag(Z)) > mtol) = nan;
            Z(Z <= 0) = nan;
            Z = real(Z);
        end

        function model = setType(model, arg, c1, c2)
            if ischar(arg)
                switch(lower(arg))
                    case {'pr', 'peng-robinson'}
                        % Peng-Robinson
                        model.eosA = 1 + sqrt(2);
                        model.eosB = 1 - sqrt(2);
                        model.omegaA = 0.4572355;
                        model.omegaB = 0.0779691;
                        model.eosType = 1;
                    case {'srk', 'soave-redlich-kwong'}
                        % Soave-Redlich-Kwong
                        model.eosA = 0;
                        model.eosB = 1;
                        model.omegaA = 0.427480;
                        model.omegaB = 0.086640;
                        model.eosType = 2;
                    case {'zj', 'zudkevitch-joffe'}
                        % Zudkevitch-Joffe
                        model.eosA = 0;
                        model.eosB = 1;
                        model.omegaA = 0.427480;
                        model.omegaB = 0.086640;
                        model.eosType = 3;
                        error('Not implemented yet.')
                    case {'rk', 'redlich-kwong'}
                        % Redlich-Kwong
                        model.eosA = 0;
                        model.eosB = 1;
                        model.omegaA = 0.427480;
                        model.omegaB = 0.086640;
                        model.eosType = 4;
                    case {'prcorr', 'peng-robinson-corrected'}
                        % Peng-Robinson
                        model.eosA = 1 + sqrt(2);
                        model.eosB = 1 - sqrt(2);
                        model.omegaA = 0.4572355;
                        model.omegaB = 0.0779691;
                        model.eosType = 5;
                    otherwise
                        error('Invalid string ''%s''.\n Valid choices are:\n PR: Peng-Robinson\n', arg);
                end
            end
            if nargin > 2
                model.eosA = c1;
                if nargin > 3
                    model.eosB = c2;
                end
            end
        end

        function Z = computeCompressibilityZ(model, p, xy, A, B, Si, Bi, isLiquid)
            if iscell(xy)
                xy = cellfun(@value, xy, 'UniformOutput', false);
                xy = [xy{:}];
            end
            if iscell(Si)
                Si = cellfun(@value, Si, 'UniformOutput', false);
                Si = [Si{:}];
            end
            if iscell(Bi)
                Bi = cellfun(@value, Bi, 'UniformOutput', false);
                Bi = [Bi{:}];
            end
            p = value(p); A = value(A); B = value(B);
                
            if model.selectGibbsMinimum
                Z0 = model.solveCubicEOS(A, B);
                bad = bsxfun(@lt, Z0, B);
                Z0(bad) = nan;

                candidates = isfinite(Z0);
                numRoots = sum(candidates, 2);
                single = numRoots == 1;
                multiple = ~single;
                Z = max(Z0, [], 2);
                if any(multiple)
                    Z_max = max(Z0(multiple, :), [], 2);
                    Z_min = min(Z0(multiple, :), [], 2);
                    xi = xy(multiple, :);

                    [~, phi] = model.computeFugacity(p(multiple), xi, Z_max, A(multiple), B(multiple), Si(multiple, :), Bi(multiple, :));
                    g_max = sum(phi.*xi, 2);
                    [~, phi] = model.computeFugacity(p(multiple), xi, Z_min, A(multiple), B(multiple), Si(multiple, :), Bi(multiple, :));
                    g_min = sum(phi.*xi, 2);
                    Zi = Z_max;
                    smallest = g_min < g_max;
                    Zi(smallest) = Z_min(smallest);
                    Z(multiple) = Zi;
                end
            else
                if isLiquid
                    Z = computeLiquidZ(model, A, B);
                else
                    Z = computeVaporZ(model, A, B);
                end
            end
        end
        
        function Z = computeLiquidZ(model, A, B)
            % Pick smallest Z factors for liquid phase (least energy)
            Z = model.solveCubicEOS(A, B);
            bad = bsxfun(@lt, Z, B);
            Z(bad) = nan;
            Z = min(Z, [], 2);
            model.checkZ(Z);
        end

        function Z = computeVaporZ(model, A, B)
            % Pick largest Z factors for vapor phase (most energy)
            Z = model.solveCubicEOS(A, B);
            bad = bsxfun(@lt, Z, B);
            Z(bad) = nan;
            Z = max(Z, [], 2);
            model.checkZ(Z);
        end
        
        function [state, report] = stepFunction(model, state, state0, dt, drivingForces, linsolve, nonlinsolve, iteration, varargin) %#ok
            % Compute a single step of the solution process for a given
            % equation of state (thermodynamic flash calculation);
            timer = tic();
            T = state.T;
            P = state.pressure;
            nc = numel(P);
            assert(all(P >= 0))
            
            z = state.components;
            ncomp = model.getNumberOfComponents();
            
            % Basic assertions
            assert(all(sum(z, 2) > 0.999), ...
                                'Molar fractions must sum up to unity.')
            assert(iteration > 0);
            [t_flash, t_stability] = deal(0);
            % If we get here, we are using a proper flash
            if iteration == 1
                % Book-keeping
                state.eos.itCount = zeros(nc, 1);
                state.eos.converged = false(nc, 1);
            end
            hasKvalues = ~isempty(model.equilibriumConstantFunctions);
            if hasKvalues
                state.K = model.evaluateEquilibriumConstants(P, T, z);
            end
            L0 = state.L;
            K0 = state.K;
            flash_method = lower(model.method);
            % Only apply calculations for cells that have not converged yet
            if iteration == 1
                x0 = state.x;
                y0 = state.y;
                % Cells which are two-phase can be flashed directly. For
                % single-phase cells, we perform a stability test in order
                % to check if the single-phase state is still stable.
                twoPhase = model.getTwoPhaseFlag(state);
                % NaN means that L was not initialized. We then start with
                % a phase stability test to get good estimates for the
                % K-values.
                twoPhase(isnan(L0)) = false;
                initSingle = ~twoPhase;
                stable = initSingle;
                % Call stability test.
                ts = tic();
                [stable(initSingle), x0(initSingle, :), y0(initSingle, :)] = ...
                    model.performPhaseStabilityTest(state.pressure(initSingle), state.T(initSingle), state.components(initSingle, :), state.K(initSingle, :));
                t_stability = toc(ts);
                updatedSingle = initSingle & ~stable;
                K0(updatedSingle, :) = y0(updatedSingle, :)./x0(updatedSingle, :);
                acf = model.CompositionalMixture.acentricFactors;
                [Si_L, Si_V, A_L, A_V, B_L, B_V, Bi] = model.getMixtureFugacityCoefficients(P, T, x0, y0, acf);
                % Solve EOS for each phase
                Z0_L = model.computeCompressibilityZ(state.pressure, x0, A_L, B_L, Si_L, Bi, true);
                Z0_V = model.computeCompressibilityZ(state.pressure, y0, A_V, B_V, Si_V, Bi, false);
                L0(~stable) = model.solveRachfordRice(L0(~stable), K0(~stable, :), z(~stable, :));
                L0(stable) = model.singlePhaseLabel(P(stable), T(stable), z(stable, :));
                active = ~stable;
                % Flag stable cells as converged
                state.eos.converged(stable) = true;
            else
                Z0_L = state.Z_L;
                Z0_V = state.Z_V;
                x0 = state.x;
                y0 = state.y;
                active = ~state.eos.converged;
            end
            if hasKvalues
                % All cells have been solved already - we are done.
                active = active | false;
                values = zeros(1, ncomp);
                resConv = true(1, ncomp);
            elseif any(active)
                state.eos.itCount(active) = state.eos.itCount(active) + 1;
                K_init = K0(active, :);
                K = K_init;
                L = L0(active);
                z = z(active, :);
                P = P(active);
                T = T(active);
                switch flash_method
                    case 'newton'
                        % Newton-based solver for equilibrium
                         updatefn = @model.newtonCompositionUpdate;
                    case 'ssi'
                        % Successive substitution solver for equilibrium
                        if iteration > model.maxSSI
                            % Remaining cells are above SSI threshold,
                            % switch to Newton
                            updatefn = @model.newtonCompositionUpdate;
                            flash_method = 'newton';
                        else
                            % Just do SSI as asked
                            updatefn = @model.substitutionCompositionUpdate;
                        end
                    otherwise
                        error('Unknown flash method %s. Try ''newton'' or ''ssi''.', model.method);
                end
                timer_flash = tic();
                [x, y, K, Z_L, Z_V, L, equilvals] = updatefn(P, T, z, K, L);
                t_flash = toc(timer_flash);
                singlePhase = L == 0 | L == 1;
                % Single phase cells are converged
                equilvals(singlePhase, :) = 0;
                values = max(equilvals, [], 1);
                conv = max(equilvals, [], 2) <= model.nonlinearTolerance;
                conv = conv & iteration > nonlinsolve.minIterations;
                resConv = values <= model.nonlinearTolerance & iteration > nonlinsolve.minIterations;
                % Insert back the local values into global arrays
                state.eos.converged(active) = conv;
                % Insert updated values in active cells
                % Just single-phase state
                K(singlePhase, :) = K_init(singlePhase, :);
                L0(active) = L;

                Z0_L(active) = value(Z_L);
                Z0_V(active) = value(Z_V);
                K0(active, :) = K;
                x0(active, :) = x;
                y0(active, :) = y;
            else
                values = zeros(1, ncomp);
                resConv = true(1, ncomp);
            end
            state.K = K0;
            state.L = L0;
            state.x = x0;
            state.y = y0;
            state.Z_L = Z0_L;
            state.Z_V = Z0_V;

            failure = false;
            failureMsg = '';
            values_converged = values <= model.nonlinearTolerance;
            if model.verbose
                printConvergenceReport(model.CompositionalMixture.names, values, values_converged, iteration);
            end
            report = model.makeStepReport(...
                            'Failure',      failure, ...
                            'FailureMsg',   failureMsg, ...
                            'Converged',    all(values_converged), ...
                            'ResidualsConverged', resConv, ...
                            'Residuals',    values);
            report.ActiveCells = sum(active);
            if model.extraOutput
                report.ActiveFlag = active;
            end
            report.Method = flash_method;
            report.TotalTime = toc(timer);
            report.FlashTime = t_flash;
            report.StabilityTime = t_stability;
        end
        
        function L = singlePhaseLabel(eos, p, T, z)
            % Li's method for phase labeling
            Vc = eos.CompositionalMixture.Vcrit;
            Tc = eos.CompositionalMixture.Tcrit;
            Vz = bsxfun(@times, Vc, z);
            
            Tc_est = sum(bsxfun(@times, Vz, Tc), 2)./sum(Vz, 2);
            L = double(T < Tc_est);
        end
        
        function [x, y, K, Z_L, Z_V, L, values] = substitutionCompositionUpdate(model, P, T, z, K, L)
            % Determine overall liquid fraction
            L = model.solveRachfordRice(L, K, z);
            % Compute liquid component fraction
            [x, sx] = model.computeLiquid(L, K, z);
            % Vapor component fraction
            y = model.computeVapor(L, K, z);
            [Z_L, Z_V, f_L, f_V] = model.getCompressibilityAndFugacity(P, T, x, y, z, [], []);
            % Compute fugacity ratios
            f_r = bsxfun(@times, sx, f_L./f_V);
            % Update equilibrium constant estimates based on fugacity ratio
            mc = model.minimumComposition;
            % Ignore tiny compositions for the purpose of updating K-values
            % and estimating convergence.
            f_r(z <= mc) = 1;
            values = abs(f_r - 1);
            K = max(K.*abs(f_r), 1e-12);
            K(~isfinite(K)) = 1;
        end

        function [stable, x, y] = performPhaseStabilityTest(model, P, T, z, K)
            if nargin < 5
                K = [];
            end
            if isempty(z)
                stable = [];
                [x, y] = deal(zeros(0, size(z, 2)));
            else
                if isempty(model.equilibriumConstantFunctions)
                    % Use fugacity-based stability test
                    z = ensureMinimumFraction(z, model.minimumComposition);
                    [stable, x, y] = phaseStabilityTest(model, z, P, T, K);
                else
                    % We have K-values, use that instead
                    K = model.evaluateEquilibriumConstants(P, T, z);
                    hasStabFn = isprop(model.PropertyModel, 'checkStabilityFunction') && ...
                                ~isempty(model.PropertyModel.checkStabilityFunction);
                    if hasStabFn
                        [stable, L] = model.PropertyModel.checkStabilityFunction(P, T, expandMatrixToCell(z));
                    else
                        L = model.solveRachfordRice([], K, z);
                        L_tol = 1e-10;
                        stable = abs(L - 1) <= L_tol | L <= L_tol;
                    end
                    x = model.computeLiquid(L, K, z);
                    y = model.computeVapor(L, K, z);
                end
            end
        end

        function K = evaluateEquilibriumConstants(model, P, T, z)
            assert(~isempty(model.equilibriumConstantFunctions), ...
                'Equilibrium constants not defined. This function is only valid when using K-value formulation.');
            isCellInput = iscell(z);
            if ~isCellInput
                z = expandMatrixToCell(z);
            end
            K = cellfun(@(fn) fn(P, T, z), model.equilibriumConstantFunctions, 'UniformOutput', false);
            
            mv = 1e8;
            for i = 1:numel(K)
                K{i}(K{i} > mv) = mv;
                K{i}(K{i} < 1/mv) = 1/mv;
            end
            if ~isCellInput
                K = value(K);
            end
        end

        
        function [x, y, K, Z_L, Z_V, L, vals] = newtonCompositionUpdate(model, P, T, z, K, L)
            [x, y, K, Z_L, Z_V, L, vals] = newtonFugacityEquilibrium(model, P, T, z, K, L);
        end

        function state = validateState(model, state)
            n_cell = size(state.pressure, 1);
            if ~isfield(state, 'L') 
                % Add empty L to indicate that we do not know. The flash
                % routine will then perform the stability test for all
                % cells.
                state.L = nan(n_cell, 1);
            end

            if ~isfield(state, 'K')
                if isfield(state, 'x') && isfield(state, 'y')
                    K = state.y./state.x;
                else
                    K = estimateEquilibriumWilson(model, state.pressure, state.T);
                end
                state.K = K;
            end
            
            if ~isfield(state, 'Z_V')
                state.Z_V = ones(n_cell, 1);
            end
            if ~isfield(state, 'Z_L')
                state.Z_L = ones(n_cell, 1);
            end
            
            if ~isfield(state, 'x')
                state.x = state.components;
            end
            if ~isfield(state, 'y')
                state.y = state.components;
            end
        end
        
        function [Pr, Tr] = getReducedPT(model, P, T, useCell)
            if nargin < 4
                useCell = true;
            end
            Tc = model.CompositionalMixture.Tcrit;
            Pc = model.CompositionalMixture.Pcrit;
            if useCell
                n = model.CompositionalMixture.getNumberOfComponents();
                Tr = cell(1, n);
                Pr = cell(1, n);
                
                Ti = 1./Tc;
                Pi = 1./Pc;
                for i = 1:n
                    Tr{i} = T.*Ti(i);
                    Pr{i} = P.*Pi(i);
                end
            else
                Tr = bsxfun(@rdivide, T, Tc);
                Pr = bsxfun(@rdivide, P, Pc);
            end
        end
        
        function [A_ij, Bi] = getMixingParameters(model, P, T, acf, useCell)
            if nargin < 5
                useCell = true;
            end
            
            % Calculate intermediate values for fugacity computation
            ncomp = model.getNumberOfComponents();
            [Pr, Tr] = model.getReducedPT(P, T, useCell);

            if useCell
                [sAi, Bi] = deal(cell(1, ncomp));
                [oA, oB] = deal(cell(1, ncomp));
                [oB{:}] = deal(model.omegaB);
            else
                oB = model.omegaB;
            end
            switch model.eosType
                case {1, 5}
                    % PR
                    if useCell
                        for i = 1:ncomp
                            ai = acf(i);
                            if model.eosType == 5 && ai > 0.49
                                oA{i} = model.omegaA.*(1 + (0.379642 + 1.48503.*ai - 0.164423.*ai.^2 + 0.016666.*ai.^3).*(1-Tr{i}.^(1/2))).^2;
                            else
                                oA{i} = model.omegaA.*(1 + (0.37464 + 1.54226.*ai - 0.26992.*ai.^2).*(1-Tr{i}.^(1/2))).^2;
                            end
                        end
                    else
                        tmp = bsxfun(@times, 0.37464 + 1.54226.*acf - 0.26992.*acf.^2, 1-Tr.^(1/2));
                        if model.eosType == 5
                            act = acf > 0.49;
                            tmp(:, act) = bsxfun(@times, 0.379642 + 1.48503.*acf(act) - 0.164423.*acf(act).^2 + 0.016666.*acf(act).^3, 1-Tr(:, act).^(1/2));
                        end
                        oA = model.omegaA.*(1 + tmp).^2;
                    end
                case 2
                    % SRK
                    if useCell
                        for i = 1:ncomp
                            oA{i} = model.omegaA.*(1 + (0.48 + 1.574.*acf(i) - 0.176.*acf(i).^2).*(1-Tr{i}.^(1/2))).^2;
                        end
                    else
                        tmp = bsxfun(@times, (0.48 + 1.574.*acf - 0.176.*acf.^2), (1-Tr.^(1/2)));
                        oA = model.omegaA.*(1 + tmp).^2;
                    end
                case 3
                    % ZJ
                    error('Not implemented yet.')
                case 4
                    % RK
                    if useCell
                        for i = 1:ncomp
                            oA{i} = model.omegaA.*Tr{i}.^(-1/2);
                        end
                    else
                        oA = model.omegaA.*Tr.^(-1/2);
                    end
                otherwise
                    error('Unknown eos type: %d', model.eosType);
            end
            bic = model.CompositionalMixture.getBinaryInteraction();
            if useCell
                A_ij = cell(ncomp, ncomp);
                for i = 1:ncomp
                    sAi{i} = ((oA{i}.*Pr{i}).^(1/2))./Tr{i};
                    Bi{i} = oB{i}.*Pr{i}./Tr{i};
                end
                for i = 1:ncomp
                    for j = i:ncomp
                        A_ij{i, j} = (sAi{i}.*sAi{j}).*(1 - bic(i, j));
                        A_ij{j, i} = A_ij{i, j};
                    end
                end
            else
                Ai = oA.*Pr./Tr.^2;
                Bi = oB.*Pr./Tr;
                A_ij = cell(ncomp, 1);
                for j = 1:ncomp
                    A_ij{j} = bsxfun(@times, bsxfun(@times, Ai, Ai(:, j)).^(1/2), 1 - bic(j, :));
                end
            end
        end
        
        function [Si_L, Si_V, A_L, A_V, B_L, B_V, Bi] = getMixtureFugacityCoefficients(model, P, T, x, y, acf)
            % Calculate intermediate values for fugacity computation
            [A_ij, Bi] = model.getMixingParameters(P, T, acf, iscell(x));
            [Si_L, A_L, B_L] = model.getPhaseMixCoefficients(x, A_ij, Bi);
            [Si_V, A_V, B_V] = model.getPhaseMixCoefficients(y, A_ij, Bi);
        end
      
        function [Z_L, Z_V, f_L, f_V] = getCompressibilityAndFugacity(model, P, T, x, y, z, Z_L, Z_V, varargin)
            [Si_L, Si_V, A_L, A_V, B_L, B_V, Bi] = model.getMixtureFugacityCoefficients(P, T, x, y, model.CompositionalMixture.acentricFactors);
            if isempty(Z_L)
                Z_L = model.computeCompressibilityZ(P, x, A_L, B_L, Si_L, Bi, true);
            end
            if isempty(Z_V)
                Z_V = model.computeCompressibilityZ(P, y, A_V, B_V, Si_V, Bi, false);
            end
            Z_L = model.setZDerivatives(Z_L, A_L, B_L, varargin{:});
            Z_V = model.setZDerivatives(Z_V, A_V, B_V, varargin{:});
            f_L = model.computeFugacity(P, x, Z_L, A_L, B_L, Si_L, Bi);
            f_V = model.computeFugacity(P, y, Z_V, A_V, B_V, Si_V, Bi);
        end
        
        function [Si, A, B] = getPhaseMixCoefficients(model, x, A_ij, Bi)
            [A, B] = deal(0);
            if iscell(x)
                ncomp = numel(x);
                Si = cell(1, ncomp);
                [Si{:}] = deal(0);
                % A = sum_ij x_i A_ij x_j
                % S_i = sum_j A_ij x_j
                % Note: A_ij is symmetric and we exploit that when we
                % evaluate the sums
                xjAij = cell(ncomp, ncomp);
                for i = 1:ncomp
                    B = B + x{i}.*Bi{i};
                    for j = 1:ncomp
                        result = A_ij{i, j}.*x{j};
                        xjAij{i, j} = result;
                        Si{i} = Si{i} + result;
                    end
                end
                for i = 1:ncomp
                    for j = i:ncomp
                        tmp = x{j}.*xjAij{j, i};
                        if i == j
                            A = A + tmp;
                        else
                            A = A + 2*tmp;
                        end
                    end
                end
            else
                ncomp = size(x, 2);
                Si = zeros(size(x));
                B = sum(x.*Bi, 2);
                A = 0;
                for i = 1:ncomp
                    A = A + sum(bsxfun(@times, bsxfun(@times, A_ij{i}, x(:, i)), x), 2);
                    Si = Si + bsxfun(@times, A_ij{i}, x(:, i));
                end
            end
        end
                
        function [E2, E1, E0] = getCubicCoefficients(model, A, B)
            [m1, m2] = model.getEOSCoefficients();

            E0 = -(A.*B + m1.*m2.*B.^2.*(B+1));
            E1  = A - (m1 + m2 - m1.*m2).*B.^2 - (m1 + m2).*B;
            E2 = (m1 + m2 - 1).*B - 1;
        end
        
        function [f, phi] = computeFugacity(model, p, x, Z, A, B, Si, Bi)
            % Compute fugacity based on EOS coefficients
            [m1, m2] = model.getEOSCoefficients();
            ncomp = model.getNumberOfComponents();
            if iscell(x)
                [f, phi] = deal(cell(1, ncomp));
                a1 = -log(Z-B);
                b1 = log((Z + m2.*B)./(Z + m1.*B)).*(A./((m1-m2).*B));
                
                Binv = 1./B;
                Ainv = 1./A;
                b2 = (Z-1).*Binv;
                for i = 1:ncomp
                    phi{i} = a1 + b1.*(2.*Si{i}.*Ainv - Bi{i}.*Binv) + Bi{i}.*b2;
                    f{i} = exp(phi{i}).*p.*x{i};
                end
            else
                a1 = -log(Z-B);
                b1 = log((Z + m2.*B)./(Z + m1.*B)).*(A./((m1-m2).*B));
                b2 = (Z-1)./B;
                tmp = bsxfun(@times, b1, 2.*bsxfun(@rdivide, Si, A) - bsxfun(@rdivide, Bi, B)) + bsxfun(@times, Bi, b2);
                phi = bsxfun(@plus, a1, tmp);
                f = exp(phi).*bsxfun(@times, p, x);
            end
        end

        function [eqs, f_L, f_V, Z_L, Z_V] = equationsEquilibrium(model, P, T, x, y, z, L, Z_L, Z_V)
            useFugacity = isempty(model.equilibriumConstantFunctions);
            ncomp = numel(x);
            if useFugacity
                [Z_L, Z_V, f_L, f_V] = model.getCompressibilityAndFugacity(P, T, x, y, z, Z_L, Z_V);
            else
                [Z_L, Z_V] = model.getCompressibilityAndFugacity(P, T, x, y, z, Z_L, Z_V);
                f_L = cell(1, ncomp);
                f_V = cell(1, ncomp);
                for i = 1:ncomp
                    K = model.equilibriumConstantFunctions{i}(P, T, z);
                    f_L{i} = K.*x{i};
                    f_V{i} = y{i};
                end
            end
            eqs = cell(1, 2*ncomp + 1);
            isLiq = value(L) == 1;
            isVap = value(L) == 0;
            isPure = isLiq | isVap;
            eqs{end} = zeros(numelValue(L), 1);
            for i = 1:ncomp
                eqs{i} = z{i} - L.*x{i} - (1-L).*y{i};
                eqs{i+ncomp} = (f_V{i} - f_L{i});
                eqs{end} = eqs{end} - ~isLiq.*y{i} + ~isVap.*x{i};
                
                if any(isPure)
                    % Treat pure phase behavior
                    eqs{i}(isPure) = z{i}(isPure) - x{i}(isPure);
                    eqs{i+ncomp}(isPure) = z{i}(isPure) - y{i}(isPure);
                    if any(isLiq)
                        eqs{end}(isLiq) = L(isLiq) - 1;
                    end
                    if any(isVap)
                        eqs{end}(isVap) = L(isVap);
                    end
                end
            end
        end

        function [x, y, L, derivatives] = getPhaseFractionAsADI(model, state, p, T, z)
            % Compute derivatives for values obtained by solving the
            % equilibrium equations (molar fractions in each phase, liquid
            % mole fraction).
            singlePhase = model.getSinglePhase(state);
            twoPhase = ~singlePhase;
            
            x = z;
            y = z;
            L = state.L;
            derivatives = struct('dL', [], 'dx', [], 'dy', []);
            if ~any(twoPhase)
                return
            end
            sample = getSampleAD(p, T, z{:});
            if isfield(state, 'FractionalDerivatives')
                fd = state.FractionalDerivatives;
                [dLdpTz, dxdpTz, dydpTz] = deal(fd.dL, fd.dx, fd.dy);
            else
                [dLdpTz, dxdpTz, dydpTz] = model.getPhaseFractionDerivativesPTZ(state, twoPhase);
            end
            ncomp = model.getNumberOfComponents();
            
            p2ph = p(twoPhase);
            T2ph = T(twoPhase);
            z2ph = cellfun(@(x) x(twoPhase), z, 'UniformOutput', false);
            
            L = model.AutoDiffBackend.convertToAD(L, sample);
            L(twoPhase) = model.addDerivativePTZ(L(twoPhase), dLdpTz, p2ph, T2ph, z2ph);
            for i = 1:ncomp
                x{i}(twoPhase) = state.x(twoPhase, i);
                y{i}(twoPhase) = state.y(twoPhase, i);
                if ~isa(z{i}, 'ADI')
                    x{i} = model.AutoDiffBackend.convertToAD(x{i}, sample);
                    y{i} = model.AutoDiffBackend.convertToAD(y{i}, sample);
                end
                x{i}(twoPhase) = model.addDerivativePTZ(x{i}(twoPhase), dxdpTz{i}, p2ph, T2ph, z2ph);
                y{i}(twoPhase) = model.addDerivativePTZ(y{i}(twoPhase), dydpTz{i}, p2ph, T2ph, z2ph);
            end
            if nargout > 3
                derivatives.dL = dLdpTz;
                derivatives.dx = dxdpTz;
                derivatives.dy = dydpTz;
            end
        end

        function [dLdpTz, dxdpTz, dydpTz] = getPhaseFractionDerivativesPTZ(model, state, cells)
            if nargin < 3
                cells = ':';
            end
            % Locally remove AD-context
            L = value(state.L);
            x = value(state.x);
            y = value(state.y);
            p = value(state.pressure);
            T = value(state.T);
            z = value(state.components);
            % Will get derivatives of these
            L = L(cells);
            x = expandMatrixToCell(x(cells, :));
            y = expandMatrixToCell(y(cells, :));
            % With respect to these
            p = p(cells);
            T = T(cells);
            z = expandMatrixToCell(z(cells, :));
            % And we need these
            Z_L = state.Z_L(cells);
            Z_V = state.Z_V(cells);
            % A few constants
            ncell = numel(p);
            ncomp = model.getNumberOfComponents();
            [xAD, yAD, zAD] = deal(cell(1, ncomp));
            
            fullZ = true;
            if fullZ
                [pAD, TAD, zAD{:}] = model.AutoDiffBackend.initVariablesAD(p, T, z{:});
                nprimary = 2 + ncomp;
            else
                [pAD, TAD, zAD{1:end-1}] = model.AutoDiffBackend.initVariablesAD(p, T, z{1:end-1});
                zAD{end} = 1;
                for i = 1:ncomp-1
                    zAD{end} = zAD{end} - zAD{i};
                end
                nprimary = 1 + ncomp;
                nprimary = nprimary + 1;
            end
            [LAD, xAD{:}, yAD{:}] = model.AutoDiffBackend.initVariablesAD(L, x{:}, y{:});

            eqsPrim  = model.equationsEquilibrium(pAD, TAD, x, y, zAD, L, Z_L, Z_V);
            eqsSec = model.equationsEquilibrium(p, T, xAD, yAD, z, LAD, Z_L, Z_V);
            ep = combineEquations(eqsPrim{:});
            es = combineEquations(eqsSec{:});
            
            % Let F ble the equilibrium equations (fugacity balance,
            % vapor/liquid balance etc). We can then write by the chain
            % rule (and using partial derivatives instead of total
            % derivatives),
            % DF / Dp = dF / dp +  dF / ds * ds/dp.
            % where p is the primary variables and s the secondary
            % variables. We then obtain 
            % ds / dp = -inv(dF / ds)*(DF / Dp)
            dFdp = ep.jac{1};
            dFds = es.jac{1};
            % We really want: dsdp = dFds\dFdp, but right hand side is
            % large so we just use LU factorization since it is much faster
            [dFds_L, dFds_U] = lu(dFds);
            dsdp = -(dFds_U\(dFds_L\dFdp));
            [I, J, V] = find(dsdp);
            % P, T and each component
            dLdpTz = zeros(ncell, nprimary);
            dxdpTz = cell(ncomp, 1);
            [dxdpTz{:}] = deal(zeros(ncell, nprimary));
            dydpTz = dxdpTz;
            
            I = I-1;
            J = J-1;
            cellNo = mod(I, ncell)+1;
            secondaryNo = floor(I/ncell)+1;
            primaryNo = floor(J/ncell)+1;
            
            isdL = secondaryNo == 1;
            makeSub = @(act) sub2ind([ncell, nprimary], cellNo(act), primaryNo(act));
            dLdpTz(makeSub(isdL)) = V(isdL);
            if ~fullZ
                dLdpTz(:, end) = -sum(dLdpTz(:, 3:end), 2);
            end
            
            % Equivialent version - if there are no zeros!
            % dLdpTz = reshape(V(isdL), ncell, nprimary);
            for i = 1:ncomp
                isdx = secondaryNo == 1 + i;
                dxdpTz{i}(makeSub(isdx)) = V(isdx);
                if ~fullZ
                    dxdpTz{i}(:, end) = -sum(dxdpTz{i}(:, 3:end), 2);
                end
                
                isdy = secondaryNo == ncomp + 1 + i;
                dydpTz{i}(makeSub(isdy)) = V(isdy);
                if ~fullZ
                    dydpTz{i}(:, end) = -sum(dydpTz{i}(:, 3:end), 2);
                end
            end
        end
        
        function L = solveRachfordRice(model, L, K, z)
            L = solveRachfordRiceVLE(L, K, z, 'min_z', model.minimumComposition);
        end

        function y = computeVapor(model, L, K, z)
            if isstruct(L)
                y = model.computeVapor(L.L, L.K, L.components);
                return
            end
            if iscell(z)
                y = z;
                sv = 0;
                for i = 1:numel(z)
                    bad = ~isfinite(K{i});
                    y{i} = K{i}.*z{i}./(L + (1-L).*K{i});
                    y{i}(bad) = 0;
                    sv = sv + value(y{i});
                end
                y = cellfun(@(k, zi) k.*zi./(L + (1-L).*k), K, z, 'UniformOutput', false);
                y = cellfun(@(x) x./sv, y, 'UniformOutput', false);
                assert(all(cellfun(@(x) all(isfinite(value(x))), y)));
            else
                y = K.*z./bsxfun(@plus, L, bsxfun(@times, 1-L, K));
                y(~isfinite(K)) = z(~isfinite(K));
                y = bsxfun(@rdivide, y, sum(y, 2));
                assert(all(isfinite(y(:))));
            end
        end
        
        function [x, sv] = computeLiquid(model, L, K, z)
            if isstruct(L)
                x = model.computeLiquid(L.L, L.K, L.components);
                return
            end
            if iscell(z)
                x = z;
                sv = 0;
                for i = 1:numel(z)
                    bad = ~isfinite(K{i});
                    x{i} = z{i}./(L + (1-L).*K{i});
                    x{i}(bad) = 0;
                    sv = sv + value(x{i});
                end
                x = cellfun(@(x) x./sv, x, 'UniformOutput', false);
                assert(all(cellfun(@(x) all(isfinite(value(x))), x)));
            else
                tmp = bsxfun(@times, 1-L, K);
                x = z./bsxfun(@plus, L, tmp);
                x(~isfinite(K)) = 0;
                sv = sum(x, 2);
                x = bsxfun(@rdivide, x, sv);
                assert(all(isfinite(x(:))));
            end
        end

        function [sL, sV] = computeSaturations(model, rhoL, rhoV, x, y, L, Z_L, Z_V)
            sL = L.*Z_L./(L.*Z_L + (1-L).*Z_V);
            sV = (1-L).*Z_V./(L.*Z_L + (1-L).*Z_V);
        end
        
        function frac = getMoleFraction(model, massfraction)
            % Convert mass fraction to molar fraction
            mm = model.CompositionalMixture.molarMass;
            if iscell(massfraction)
                ncomp = numel(massfraction);
                moles = cell(1, ncomp);
                totMole = 0;
                for i = 1:numel(massfraction)
                    moles{i} = massfraction{i}./mm(i);
                    totMole = totMole + moles{i};
                end
                frac = cellfun(@(x) x./totMole, moles, 'UniformOutput', false);
            else
                moles = bsxfun(@times, massfraction, 1./mm);
                frac = bsxfun(@rdivide, moles, sum(moles, 2));
            end
        end
        
        function frac = getMassFraction(model, molfraction)
            % Convert molar fraction to mass fraction
            if iscell(molfraction)
                ncomp = numel(molfraction);
                mass = cell(1, ncomp);
                totMass = 0;
                for i = 1:ncomp
                    mi = molfraction{i};
                    if ~isempty(mi)
                        mass{i} = model.CompositionalMixture.molarMass(i).*molfraction{i};
                        totMass = totMass + mass{i};
                    end
                end
                frac = cell(size(mass));
                for i = 1:ncomp
                    if ~isempty(mass{i})
                        frac{i} = mass{i}./totMass;
                    end
                end
            else
                mass = bsxfun(@times, molfraction, model.CompositionalMixture.molarMass);
                frac = bsxfun(@rdivide, mass, sum(mass, 2));
            end
        end

        function L = estimateSinglePhaseState(model, p, T, z, L, stable)
            z = z(stable, :);
            p = p(stable);
            T = T(stable);

            K = estimateEquilibriumWilson(model, p, T);
            L0 = repmat(0.5, nnz(stable), 1);
            L(stable) = model.solveRachfordRice(L0, K, z);

            L(stable & L > 0.5) = 1;
            L(stable & L <= 0.5) = 0;
        end

        function Z = setZDerivatives(model, Z, A, B, cellJacMap)
            % Z comes from the solution of a cubic equation of state, so
            % the derivatives are not automatically computed. By
            % differentiating the cubic EOS manually and solving for dZ/du
            % where u is some primary variable, we can still obtain
            % derivatives without making any assumptions other than the EOS
            % being a cubic polynomial
            S = getSampleAD(A, B);
            if isnumeric(S)
                % This function does nothing if A, B are both doubles.
                return;
            end

            backend = model.AutoDiffBackend;
            isDiagonal = isa(backend, 'DiagonalAutoDiffBackend');
            if isnumeric(Z)
                Z = backend.convertToAD(Z, S);
            end
            if nargin < 5
                if isDiagonal
                    cellJacMap = cell(Z.offsets(end)-1, 1);
                else
                    cellJacMap = cell(numel(Z.jac), 1);
                end
            end
            [E2, E1, E0] = model.getCubicCoefficients(A, B);
            e2 = value(E2);
            e1 = value(E1);
            z = value(Z);
            if isDiagonal
                offset = Z.offsets;
                rowMajor = backend.rowMajor;
                if numel(Z.jac)
                    % Should be diagonal
                    allow_implicit = Z.jac{1}.allowImplicitExpansion;
                else
                    allow_implicit = false; % Does not matter.
                end
                if allow_implicit
                    % Newer Matlab allows implicit expansion - use that, if
                    % available
                    mlt = @(x, y) x.*y;
                    div = @(x, y) x./y;
                else
                    mlt = @(x, y) bsxfun(@times, x, y);
                    div = @(x, y) bsxfun(@rdivide, x, y);
                end
                if rowMajor
                    z = z';
                    e1 = e1';
                    e2 = e2';
                end
                for i = 1:numel(Z.jac)
                    act = offset(i):offset(i+1)-1;
                    map = cellJacMap(act);
                    map = map(~cellfun(@isempty, map));

                    if isempty(map)
                        if Z.jac{i}.dim(1) ~= numel(z) && numel(Z.jac{i}.subset) ~= numel(z)
                            continue
                        end
                        % All cells
                        dE2 = E2.jac{i}.diagonal;
                        dE1 = E1.jac{i}.diagonal;
                        dE0 = E0.jac{i}.diagonal;

                        % d = -(dE2.*z.^2 + dE1.*z + dE0)./(3*z.^2 + 2*z.*e2 + e1);
                        d = -div(mlt(dE2, z.^2) + mlt(dE1, z) + dE0, 3*z.^2 + 2*z.*e2 + e1);
                        if any(any(d~=0))
                            if Z.jac{i}.isZero
                                Z.jac{i} = Z.jac{i}.expandZero();
                            end
                            Z.jac{i}.diagonal = d;
                            if all(Z.jac{i}.subset == 0)
                                Z.jac{i}.subset = [];
                            end
                        end
                    else
                        % Subset of cells
                        map = map{1};
                        dE2 = E2.jac{i}.getDiagonal(map);
                        dE1 = E1.jac{i}.getDiagonal(map);
                        dE0 = E0.jac{i}.getDiagonal(map);
                        if isempty(dE0)
                            assert(isempty(dE1));
                            assert(isempty(dE2));
                            continue
                        end
                        zi = z(map);
                        e2i = e2(map);
                        e1i = e1(map);
                        d = -(dE2.*zi.^2 + dE1.*zi + dE0)./(3*zi.^2 + 2*zi.*e2i + e1i);
                        if any(any(d~=0))
                            if Z.jac{i}.isZero
                                Z.jac{i} = Z.jac{i}.expandZero();
                            end
                            Z.jac{i} = Z.jac{i}.setDiagonal(map, d);
                            if ~isempty(Z.jac{i}.subset)
                                Z.jac{i}.subset(map) = (1:numel(map))';
                            end
                        end
                    end
                end
            else
                for i = 1:numel(Z.jac)
                    if isempty(cellJacMap{i})
                        [n, m] = size(Z.jac{i});
                        if n == 0 || m ~= numel(z)
                            % There are either no derivatives present or the
                            % derivatives are not of the right dimension
                            continue
                        end
                        dE2 = model.getJac(E2, i);
                        dE1 = model.getJac(E1, i);
                        dE0 = model.getJac(E0, i);

                        d = -(dE2.*z.^2 + dE1.*z + dE0)./(3*z.^2 + 2*z.*e2 + e1);
                        Z.jac{i} = sparse((1:n)', (1:m)', d, n, m);
                    else
                        map = cellJacMap{i};
                        zi = z(map);
                        e2i = e2(map);
                        e1i = e1(map);
                        
                        [n, m] = size(Z.jac{i});
                        if n == 0 || m ~= numel(zi)
                            % There are either no derivatives present or the
                            % derivatives are not of the right dimension
                            continue
                        end
                        dE2 = model.getJacSub(E2, i, map);
                        dE1 = model.getJacSub(E1, i, map);
                        dE0 = model.getJacSub(E0, i, map);

                        d = -(dE2.*zi.^2 + dE1.*zi + dE0)./(3*zi.^2 + 2*zi.*e2i + e1i);
                        Z.jac{i} = sparse(map, (1:m)', d, n, m);
                    end
                end
            end
        end

        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@PhysicalModel(model, state0, state, dt, drivingForces);
            changed = false;
            state = model.setFlag(state);
            if ~isfield(state, 'x')
                state.x = model.computeLiquid(state);
                changed = true;
            end
            if ~isfield(state, 'y')
                state.y = model.computeVapor(state);
                changed = true;
            end
            if changed
                singlePhase = state.L == 1 | state.L == 0;
                z_1ph = state.components(singlePhase, :);
                state.x(singlePhase, :) = z_1ph;
                state.y(singlePhase, :) = z_1ph;
            end
        end

        function state = setFlag(model, state, pureLiquid, pureVapor)
            if nargin < 4
                L = value(state.L);
                pureVapor = L == 0;
                if nargin < 3
                    pureLiquid = L == 1;
                end
            end

            pureVapor(pureLiquid) = false;
            state.flag = 1*pureLiquid + 2*pureVapor;
        end

        function isSinglePhase = getSinglePhase(model, state)
            if ~isfield(state, 'flag')
                state = model.setFlag(state);
            end
            isSinglePhase = state.flag ~= 0;
        end

        function [isLiquid, isVapor, is2ph] = getFlag(model, state)
            if ~isfield(state, 'flag')
                state = model.setFlag(state);
            end
            flag = state.flag;
            isLiquid = flag == 1;
            isVapor = flag == 2;
            is2ph = ~(isLiquid | isVapor);
        end
        
        function is2ph = getTwoPhaseFlag(model, state)
            if ~isfield(state, 'flag')
                state = model.setFlag(state);
            end
            is2ph = state.flag == 0;
        end
    end

    methods(Access = protected)
        function [m1, m2] = getEOSCoefficients(model)
            m1 = model.eosA;
            m2 = model.eosB;
        end
    end
    
    methods (Static, Access=protected)
        function dx = getJac(x, ix)
            if isa(x, 'ADI')
                dx = diag(x.jac{ix});
            else
                dx = 0;
            end
        end

        function dx = getJacSub(x, ix, map)
            if isa(x, 'ADI')
                [ii, jj, vv] = find(x.jac{ix});
                if isempty(vv)
                    dx = 0;
                else
                    keep = ii == map(jj);
                    dx = zeros(size(map));
                    dx(jj(keep)) = vv(keep);
                end
            else
                dx = 0;
            end
        end

        function D = makeDiagonal(x, vx, ncell)
            D = sparse(vx, vx, x, ncell, ncell);
        end

        function checkZ(Z)
            % Throw error if Z takes on unphysical values: Negative, or
            % non-finite values.
            if any(Z < 0)
                error('%d negative Z-factors detected...', sum(Z < 0));
            end
            if any(~isfinite(Z))
                error('%d non-finite Z-factors detected...', sum(~isfinite(Z)));
            end
        end

        function Z = fallbackCubic(A, E2, E1, E0)
            n = numel(A);
            Z = zeros(n, 3);
            for i = 1:n
                Z(i, :) = roots([1, E2(i), E1(i), E0(i)]);
            end
        end
        
        function result = addDerivativePTZ(result, dvdptz, p, T, z)
            if ~isa(result, 'ADI')
                return
            end
            nJac = numel(result.jac);
            for jacNo = 1:nJac
                result.jac{jacNo} = 0*result.jac{jacNo};
            end
            variables = [{p}, {T}, z];
            for varNo = 1:numel(variables)
                v = variables{varNo};
                if isa(v, 'ADI')
                    v_der = dvdptz(:, varNo);
                    D = [];
                    for jacNo = 1:nJac
                        [inc, D] = diagMult(v_der, v.jac{jacNo}, D);
                        result.jac{jacNo} = result.jac{jacNo} + inc;
                    end
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


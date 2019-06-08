classdef EquationOfStateModel < PhysicalModel
    % Equation of state model. Implements generalized two-parameter cubic
    % equation of state with Newton and successive substitution solvers, as
    % well as standard functions for computing density and viscosity.
    properties
        fluid % CompositionalFluid
        omegaA % Parameter for EOS
        omegaB % Parameter for EOS
        useNewton % Use Newton based solver for flash. If set to false, successive substitution if used instead.
        PropertyModel % Model to be used for property evaluations
        selectGibbsMinimum = true; % Use minimum Gibbs energy to select Z
        alpha = [];
        minimumComposition = 1e-8; % Minimum composition value (for numerical stability)
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
            assert(isa(fluid, 'CompositionalFluid'));
            model.fluid = fluid;
            model.nonlinearTolerance = 1e-4;
            model.useNewton = false;
            model.fastDerivatives = true;
            model.PropertyModel = CompositionalPropertyModel(fluid);
            model = model.setType(eosname);
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

        function [m1, m2] = getEOSCoefficients(model)
            m1 = model.eosA;
            m2 = model.eosB;
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
                
            if ~model.selectGibbsMinimum
                if isLiquid
                    Z = computeLiquidZ(model, A, B);
                else
                    Z = computeVaporZ(model, A, B);
                end
                return
            end
            
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
            T = state.T;
            P = state.pressure;
            nc = numel(P);
            assert(all(P >= 0))
            
            z = state.components;
            ncomp = model.fluid.getNumberOfComponents();
            
            % Basic assertions
            assert(all(sum(z, 2) > 0.999), ...
                                'Molar fractions must sum up to unity.')
            assert(iteration > 0);
            
            if iteration == 1
                % Book-keeping
                state.eos.itCount = zeros(nc, 1);
                state.eos.converged = false(nc, 1);
            end

            L0 = state.L;
            K0 = state.K;
            % Only apply calculations for cells that have not converged yet
            if iteration == 1
                x0 = state.x;
                y0 = state.y;
                [~, ~, twoPhase] = model.getFlag(state);
                initSingle = ~twoPhase;
                stable = initSingle;
                [stable(initSingle), x0(initSingle, :), y0(initSingle, :)] = model.performPhaseStabilityTest(state.pressure(initSingle), state.T(initSingle), state.components(initSingle, :));
                acf = model.fluid.acentricFactors;
                [Si_L, Si_V, A_L, A_V, B_L, B_V, Bi] = model.getMixtureFugacityCoefficients(P, T, x0, y0, acf);
                % Solve EOS for each phase
                Z0_L = model.computeCompressibilityZ(state.pressure, x0, A_L, B_L, Si_L, Bi, true);
                Z0_V = model.computeCompressibilityZ(state.pressure, y0, A_V, B_V, Si_V, Bi, false);
                L0 = model.solveRachfordRice(L0, K0, z);
                L0(stable & L0 >  0.5) = 1;
                L0(stable & L0 <= 0.5) = 0;
                % L0 = model.estimateSinglePhaseState(state.pressure, state.T, state.components, L0, stable);
                active = ~stable;
            else
                Z0_L = state.Z_L;
                Z0_V = state.Z_V;
                x0 = state.x;
                y0 = state.y;
                active = ~state.eos.converged;
            end
            
            if any(active)
                state.eos.itCount(active) = state.eos.itCount(active) + 1;
                L = L0(active);
                z = z(active, :);
                K = K0(active, :);
                P = P(active);
                T = T(active);

                if model.useNewton
                    % Newton-based solver for equilibrium
                    [x, y, K, Z_L, Z_V, L, equilvals] = model.newtonCompositionUpdate(P, T, z, K, L);
                else
                    % Successive substitution solver for equilibrium
                    [x, y, K, Z_L, Z_V, L, equilvals] = model.substitutionCompositionUpdate(P, T, z, K, L);
                end

                values = max(equilvals, [], 1);
                valconv = values <= model.nonlinearTolerance;
                conv = max(equilvals, [], 2) <= model.nonlinearTolerance;
                conv = conv & iteration > nonlinsolve.minIterations;
                resConv = values <= model.nonlinearTolerance & iteration > nonlinsolve.minIterations;
                % Insert back the local values into global arrays
                state.eos.converged(active) = conv;
                % Insert updated values in active cells
                L0(active) = L;

                Z0_L(active) = Z_L;
                Z0_V(active) = Z_V;
                K0(active, :) = K;
                x0(active, :) = x;
                y0(active, :) = y;
            else
                valconv = true(1, ncomp);
                values = zeros(1, ncomp);
                resConv = true(1, ncomp);
            end
            state.K = K0;
            state.L = L0;
            state.x = x0;
            state.y = y0;
            state.Z_L = Z0_L;
            state.Z_V = Z0_V;

            state.mixing.K = state.K;

            failure = false;
            failureMsg = '';
            
            if model.verbose
                printConvergenceReport(model.fluid.names, values, valconv, iteration);
            end
            report = model.makeStepReport(...
                            'Failure',      failure, ...
                            'FailureMsg',   failureMsg, ...
                            'Converged',    all(valconv), ...
                            'ResidualsConverged', resConv, ...
                            'Residuals',    values);
            report.ActiveCells = sum(active);
        end
        
        function [x, y, K, Z_L, Z_V, L, values] = substitutionCompositionUpdate(model, P, T, z, K, L)
            % Determine overall liquid fraction
            L = model.solveRachfordRice(L, K, z);
            % Compute liquid component fraction
            [x, sx] = model.computeLiquid(L, K, z);
            % Vapor component fraction
            y = model.computeVapor(L, K, z);
            [Z_L, Z_V, f_L, f_V] = model.getCompressibilityAndFugacity(P, T, x, y, z, [], []);
            
            % Pure liquid / pure vapor is converged
            ok = L == 1 | L == 0;
            % Compute fugacity ratios
            f_r = bsxfun(@times, sx, f_L./f_V);
            f_r(ok, :) = 1;
            f_r(z == 0) = 1;
            % Update equilibrium constant estimates based on fugacity ratio
            values = abs(f_r - 1);
            K = max(K.*abs(f_r), 1e-12);
            K(~isfinite(K)) = 1;
        end

        function [stable, x, y] = performPhaseStabilityTest(model, P, T, z, cells)
            if nargin < 5
                cells = [];
            end
            if isempty(z)
                stable = [];
                [x, y] = deal(zeros(0, size(z, 2)));
            else
                z = ensureMinimumFraction(z, model.minimumComposition);
                [stable, x, y] = phaseStabilityTest(model, z, P, T, z, z);
            end
        end

        function [x, y, K, Z_L, Z_V, L, vals] = newtonCompositionUpdate(model, P, T, z, K, L)
            [x, y, K, Z_L, Z_V, L, vals] = newtonFugacityEquilibrium(model, P, T, z, K, L);
        end

        function state = validateState(model, state)
            n_cell = size(state.pressure, 1);
            if ~isfield(state, 'L') 
                % Initial guess of 0.5
                state.L = repmat(0.5, n_cell, 1);
            end

            if ~isfield(state, 'K')
                state.K = estimateEquilibriumWilson(model, state.pressure, state.T);
            else
                pure = state.L == 1 | state.L == 0;
                if any(pure)
                    state.K(pure, :) = estimateEquilibriumWilson(model, state.pressure(pure), state.T(pure));
                end
            end
            
            if ~isfield(state, 'Z_V')
                state.Z_V = ones(n_cell, 1);
            end
            if ~isfield(state, 'Z_V')
                state.Z_V = ones(n_cell, 1);
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
            
            if useCell
                n = model.fluid.getNumberOfComponents();
                Tr = cell(1, n);
                Pr = cell(1, n);

                for i = 1:n
                    Tr{i} = T./model.fluid.Tcrit(i);
                    Pr{i} = P./model.fluid.Pcrit(i);
                end
            else
                Tr = bsxfun(@rdivide, T, model.fluid.Tcrit);
                Pr = bsxfun(@rdivide, P, model.fluid.Pcrit);
            end
        end
        
        function [A_ij, Bi] = getMixingParameters(model, P, T, acf, useCell)
            if nargin < 5
                useCell = true;
            end
            
            % Calculate intermediate values for fugacity computation
            ncomp = model.fluid.getNumberOfComponents();
            [Pr, Tr] = model.getReducedPT(P, T, useCell);

            if useCell
                [Ai, Bi] = deal(cell(1, ncomp));
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
                            if model.eosType == 5 && acf(i) > 0.49
                                oA{i} = model.omegaA.*(1 + (0.379642 + 1.48503.*acf(i) - 0.164423.*acf(i).^2 + 0.016666.*acf(i).^3).*(1-Tr{i}.^(1/2))).^2;
                            else
                                oA{i} = model.omegaA.*(1 + (0.37464 + 1.54226.*acf(i) - 0.26992.*acf(i).^2).*(1-Tr{i}.^(1/2))).^2;
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
            bic = model.fluid.getBinaryInteraction();
            if useCell
                A_ij = cell(ncomp, ncomp);
                for i = 1:ncomp
                    tmp = Pr{i}./Tr{i};
                    Ai{i} = oA{i}.*tmp./Tr{i};
                    Bi{i} = oB{i}.*tmp;
                end
                for i = 1:ncomp
                    for j = i:ncomp
                        A_ij{i, j} = (Ai{i}.*Ai{j}).^(1/2).*(1 - bic(i, j));
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
      
        function [Z_L, Z_V, f_L, f_V] = getCompressibilityAndFugacity(model, P, T, x, y, z, sO, sG, state, varargin)
            if nargin < 9
                state = [];
            end
            [Si_L, Si_V, A_L, A_V, B_L, B_V, Bi] = model.getMixtureFugacityCoefficients(P, T, x, y, model.fluid.acentricFactors);
            if isfield(state, 'Z_L')
                Z_L = state.Z_L;
            else
                Z_L = model.computeCompressibilityZ(P, x, A_L, B_L, Si_L, Bi, true);
            end
            if isfield(state, 'Z_V')
                Z_V = state.Z_V;
            else
                Z_V = model.computeCompressibilityZ(P, y, A_V, B_V, Si_V, Bi, false);
            end
            if iscell(x)
                s = getSampleAD(P, T, x{:}, y{:});
                if isa(s, 'GenericAD')
                    Z_L = double2GenericAD(Z_L, s);
                    Z_V = double2GenericAD(Z_V, s);
                elseif isa(s, 'ADI')
                    Z_L = double2ADI(Z_L, s);
                    Z_V = double2ADI(Z_V, s);
                end
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
                for i = 1:ncomp
                    B = B + x{i}.*Bi{i};
                    for j = 1:ncomp
                        A_ijx_j = A_ij{i, j}.*x{j};
                        A = A + x{i}.*A_ijx_j;
                        Si{i} = Si{i} + A_ijx_j;
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
            ncomp = model.fluid.getNumberOfComponents();
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
            t1 = tic();
            [Z_L, Z_V, f_L, f_V] = model.getCompressibilityAndFugacity(P, T, x, y, z, Z_L, Z_V);
            
            ncomp = numel(x);
            timer = tic();
            eqs = cell(1, 2*ncomp + 1);
            sample = getSampleAD(P, T, x{:}, y{:}, z{:});
            emptyJac = double2ADI(zeros(size(value(P))), sample);
            
            isLiq = value(L) == 1;
            isVap = value(L) == 0;
            isPure = isLiq | isVap;
            eqs{end} = emptyJac + (isLiq | isVap);
            for i = 1:ncomp
                eqs{i} = z{i} - L.*x{i} - (1-L).*y{i} + emptyJac;
                eqs{i+ncomp} = (f_V{i} - f_L{i}) + emptyJac;
                eqs{end} = eqs{end} - ~isLiq.*y{i} + ~isVap.*x{i} + emptyJac;
                
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

        function [x, y, LL] = getPhaseFractionAsADI(model, state, pP, TP, z0)
            % Compute derivatives for values obtained by solving the
            % equilibrium equations (molar fractions in each phase, liquid
            % mole fraction).
            twoPhase = true(size(value(pP)));

            % This function currently only works for the full set of
            % variables, making the local systems a bit too large. The root
            % cause of this is that the derivatives must be massaged to
            % take subsets.
            
            % [~, ~, twoPhase] = model.getFlag(state);
            % twoPhase = state.L > 0 & state.L < 1;

            [x, y] = deal(z0);
            LL = state.L;
            
            primVar = getSampleAD(pP, TP, z0{:});
            for i = 1:numel(x)
                if ~isa(x{i}, 'ADI')
                    x{i} = double2ADI(x{i}, primVar);
                    y{i} = double2ADI(y{i}, primVar);
                end
            end
            if ~isa(LL, 'ADI')
                LL = double2ADI(LL, primVar);
            end
            if ~any(twoPhase)
                return
            end

            z0 = cellfun(@(x) x(twoPhase), z0, 'UniformOutput', false);
            pP = pP(twoPhase);
            TP = TP(twoPhase);
            zP = z0;
            
            L = state.L(twoPhase);
            K = state.K(twoPhase, :);
            zD = cellfun(@value, zP, 'UniformOutput', false);
            % Compute secondary variables as pure doubles
            xD = expandMatrixToCell(state.x(twoPhase, :));
            yD = expandMatrixToCell(state.y(twoPhase, :));
            ncomp = numel(xD);
            Z_L = state.Z_L(twoPhase);
            Z_V = state.Z_V(twoPhase);

            xS = cell(1, ncomp);
            yS = cell(1, ncomp);
            zS = zD;
            
            % Secondary as primary variables to get equilibrium jacobian
            % with respect to them.
            [xS{:}, yS{:}, LS] = initVariablesADI(xD{:}, yD{:}, L);

            xP = xD;
            yP = yD;

            pS = value(pP);
            TS = value(TP);
            LP = L;
            
            eqsPrim  = model.equationsEquilibrium(pP, TP, xP, yP, zP, LP, Z_L, Z_V);
            eqsSec = model.equationsEquilibrium(pS, TS, xS, yS, zS, LS, Z_L, Z_V);
            
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

            % We now have a big lumped matrix of all the derivatives, which
            % we can then extract into the regular ADI format of cell
            % arrays
            [I, J, V] = find(dsdp);
            jSz = cellfun(@(x) size(x, 2), primVar.jac);
            V(~isfinite(V)) = 0;

            jOffset = cumsum([0, jSz]);
            jump = numel(L);
            nj = numel(primVar.jac);
            for i = 1:ncomp
                startIx = jump*(i-1);
                endIx = jump*i;
                startIy = startIx + ncomp*jump;
                endIy = endIx + ncomp*jump;
                for j = 1:nj
                    startJ = jOffset(j);
                    endJ = jOffset(j+1);

                    actx = I > startIx & I <= endIx & J > startJ & J <= endJ;
                    acty = I > startIy & I <= endIy & J > startJ & J <= endJ;
                    
                    x{i}.jac{j}(twoPhase, :) = sparse(I(actx) - startIx, J(actx) - startJ, V(actx), endIx-startIx, endJ-startJ);
                    y{i}.jac{j}(twoPhase, :) = sparse(I(acty) - startIy, J(acty) - startJ, V(acty), endIy-startIy, endJ-startJ);
                end
                x{i}.val(twoPhase) = xD{i};
                y{i}.val(twoPhase) = yD{i};
            end
            startI = 2*ncomp*jump;
            endI = startI + jump;
            for i = 1:nj
                startJ = jOffset(i);
                endJ = jOffset(i+1);
                actL = I > startI & J > startJ & J <= endJ;
                LL.jac{i}(twoPhase, :) = sparse(I(actL) - startI, J(actL) - startJ, V(actL), endI-startI, endJ-startJ);
            end
        end

        function L = solveRachfordRice(model, L, K, z) %#ok
            L = solveRachfordRiceVLE(L, K, z);
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
            if iscell(massfraction)
                ncomp = numel(massfraction);
                moles = cell(1, ncomp);
                totMole = 0;
                for i = 1:numel(massfraction)
                    moles{i} = massfraction{i}./model.fluid.molarMass(i);
                    totMole = totMole + moles{i};
                end
                frac = cellfun(@(x) x./totMole, moles, 'UniformOutput', false);
            else
                moles = bsxfun(@times, massfraction, 1./model.fluid.molarMass);
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
                        mass{i} = model.fluid.molarMass(i).*molfraction{i};
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
                mass = bsxfun(@times, molfraction, model.fluid.molarMass);
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
            if nargin < 5
                if isa(Z, 'GenericAD')
                    cellJacMap = cell(Z.offsets(end)-1, 1);
                elseif isa(Z, 'ADI')
                    cellJacMap = cell(numel(Z.jac), 1);
                else
                    cellJacMap = {};
                end
            end
            [E2, E1, E0] = model.getCubicCoefficients(A, B);
            e2 = value(E2);
            e1 = value(E1);
            z = value(Z);
            if isa(Z, 'GenericAD')
                offset = Z.offsets;
                for i = 1:numel(Z.jac)
                    act = offset(i):offset(i+1)-1;
                    map = cellJacMap(act);
                    map = map(~cellfun(@isempty, map));

                    if isempty(map)
                        if Z.jac{i}.dim(1) ~= numel(z)
                            continue
                        end
                        % All cells
                        dE2 = E2.jac{i}.diagonal;
                        dE1 = E1.jac{i}.diagonal;
                        dE0 = E0.jac{i}.diagonal;

                        d = -(dE2.*z.^2 + dE1.*z + dE0)./(3*z.^2 + 2*z.*e2 + e1);
                        if any(any(d~=0))
                            if Z.jac{i}.isZero
                                Z.jac{i} = Z.jac{i}.expandZero();
                            end
                            Z.jac{i}.diagonal = d;
                            Z.jac{i}.subset = [];
                        end
                    else
                        % Subset of cells
                        map = map{1};
                        dE2 = E2.jac{i}.diagonal(map, :);
                        dE1 = E1.jac{i}.diagonal(map, :);
                        dE0 = E0.jac{i}.diagonal(map, :);
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
                            Z.jac{i}.diagonal(map, :) = d;
                            if ~isempty(Z.jac{i}.subset)
                                Z.jac{i}.subset(map) = (1:numel(map))';
                            end
                        end
                    end
                end
                
            elseif isa(Z, 'ADI')
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
            state = model.setFlag(state);
            state.x = model.computeLiquid(state);
            state.y = model.computeVapor(state);
        end

        function state = setFlag(model, state, pureLiquid, pureVapor)
            if nargin < 4
                pureVapor = state.L == 0;
                if nargin < 3
                    pureLiquid = state.L == 1;
                end
            end
            if size(state, 2) > 2
                pureWater = sum(state.s(:, 2:end), 2) == 0;
                pureVapor(pureWater) = false;
                pureLiquid(pureWater) = true;
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
    end

    methods (Static, Access=private)
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
    end
end

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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


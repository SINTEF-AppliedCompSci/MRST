classdef NumericalPressureReductionFactors < StateFunction
    properties
        useDiagonalReduction
        includeDerivatives = true;
        fullDerivatives    = false;
        checkDerivatives   = false;
        rowMajor           = false;
    end
    methods
        
        function prf = NumericalPressureReductionFactors(model, varargin)
            prf@StateFunction(model, varargin{:});
            AD = model.AutoDiffBackend;
            isDiag = isa(AD, 'DiagonalAutoDiffBackend');
            prf.useDiagonalReduction = isDiag && AD.useMex;
            if isDiag
                prf.rowMajor = AD.rowMajor;
            end
        end
        
        function weights = evaluateOnDomain(prop, model, state)
            % Get state with the full primary variable set as AD
            [stateAD, ncell, ncomp] = prop.getStateAD(model, state);
            % Timestep and mass at current and previous timestep
            dt    = state.reductionFactorProps.dt;
            % Pressure at current and previous iteration
            p      = value(state.pressure);
            mass  = model.getProps(stateAD, 'ComponentTotalMass');            
            factor = state.reductionFactorProps;
            factor.mass = mass;
            factor.pressure0 = factor.pressure;
            factor.pressure = p;
            % Weights at current and previous iteration
            w      = prop.getWeights(factor, dt, ncell, ncomp);
            w_prev = state.reductionFactorProps.weights;
            % Get derivatives of weights wrt pressure
            % Construct AD weights
            pAD = state.pressure;
            isAD = isa(pAD, 'ADI') && prop.includeDerivatives;
            if isAD
                dwdp = prop.getWeightDerivatives(w, w_prev, factor);
                hasWeights = ~isempty(dwdp);
            else
                hasWeights = false;
            end
            weights = cell(ncomp, 1);
            for i = 1:ncomp
                Wp = w(:, i);
                if hasWeights
                    dpi = dwdp(:, i);
                    if any(dpi)
                        Wp = model.AutoDiffBackend.convertToAD(Wp, pAD);
                        Wp.jac{1} = sparse(1:ncell, 1:ncell, dpi, ncell, ncell);
                    end
                end
                weights{i} = Wp;
            end
            if prop.checkDerivatives && isAD
                tmp = 0;
                for i = 1:ncomp
                    m = mass{i};
                    Wp = w(:, i);
                    if hasWeights
                        dpi = dwdp(:, i);
                        if any(dpi)
                            Wp = model.AutoDiffBackend.convertToAD(Wp, m);
                            if isa(Wp.jac{1}, 'DiagonalJacobian')
                                Wp.jac{1}.diagonal(:, 1) = dpi;
                            else
                                sz = ncell;
                                Wp.jac{1} = sparse(1:ncell, 1:ncell, dpi, ncell, sz);
                            end
                            
                        end
                    end
                    tmp = tmp + m.*Wp;
                end
                for derNo = 2:ncomp
                    d = getDiagonal(tmp, derNo);
                    bad = abs(d) > 1e-8*mean(abs(getDiagonal(tmp, 1)));
                    if any(bad)
                        warning('A_p with respect to variable %d has %d non-zero derivatives.', derNo, sum(bad));
                    end
                end
            end
        end
        
        function w = getWeights(prp, factor, dt, ncell, ncomp)
            % Compute mass difference
            mass = factor.mass;
            mass0 = factor.mass0;
            acc = cellfun(@(m, m0) (m - m0)./dt, mass, mass0, 'UniformOutput', false);
            isD = cellfun(@isnumeric, acc);
            if all(isD)
                w = ones(ncell, ncomp);
                return;
            end
            ndof = ncell*ncomp;
            if prp.useDiagonalReduction
                % Special case: We can assume that all jacobians are diagonal.
                diags = cellfun(@(x) x.jac{1}.diagonal, acc, 'UniformOutput', false);
                extra = cell(1, ncomp);
                b = zeros(ndof, 1);
                scale = 1/barsa;
                % scale = sum(M(1:ncomp:end-ncomp+1, :), 1);
                b(1:ncomp:end-ncomp+1) = scale;
                if prp.fullDerivatives && ~isempty(factor.weights)
                    masses = value(mass');
                    masses0 = value(mass0');
                    dm = masses - masses0;
                    dm(abs(dm) < 1e-8) = 0;
                    % Matrix entries
                    maindiag = zeros(ncell, ncomp);
                    for i = 1:ncomp
                        d = getDiagonal(mass{i}, i);
                        e = masses.*d./dm(:, i);
                        e(~isfinite(e)) = 0;
                        extra{i} = e;
                        maindiag(:, i) = d;
                    end
                    % Right-hand side
                    w0 = sum(factor.weights, 2);
                    b_extra = bsxfun(@times, maindiag./dm, w0);
                    b_extra(~isfinite(b_extra)) = 0;
                    b_extra = b_extra';
                    b = b + b_extra(:);
                    for i = 1:ncomp
                        diags{i} = diags{i} + extra{i};
                    end
                end
                diags = cellfun(@(x) x, diags, 'UniformOutput', false);
                if ~prp.rowMajor
                    diags = cellfun(@(x) x', diags, 'UniformOutput', false);
                end
                
                M = vertcat(diags{:});
                ncell = size(M, 2);
                sz = repmat(ncomp, ncell, 1);
                Mi = invv(reshape(M, [], 1), sz);
                [i, j] = blockDiagIndex(sz, sz);
                M_inv = sparse(i, j, Mi, ndof, ndof);
                w = M_inv*b;
                w = reshape(w, ncomp, [])';
            else
                c = combineEquations(acc{:});
                % Get mass difference Jacobians
                A = c.jac{1};
                % Safeguard against singular system (typically incomp/ weakly incomp)
                [~, jj, v] = find(A);
                ix  = jj <= ncell;
                tol = 1e-16;
                if ~any(ix) || norm(v(ix), inf) < tol
                    A = A + sparse(1:ncell, 1:ncell, 1, ndof, ndof);
                end
                b = zeros(1, ndof);
                b(1:ncell) = 1/barsa;
                % Compute weights
                w = b/A;
                w = reshape(w', [], ncomp);
            end
        end
        
        function [state, ncell, ncomp] = getStateAD(prop, model, state)
            % Remove existing AD
            state = model.reduceState(state, true);
            % Get variables
            [vars, names, origin] = model.getPrimaryVariables(state);
            isP   = strcmp(names, 'pressure');
            origP = origin{isP};
            isW   = cellfun(@(x) ~strcmp(x, origP), origin);
            % Initialize AD
            [vars{~isW}] = model.AutoDiffBackend.initVariablesAD(vars{~isW});
            state = model.initStateAD(state, vars, names, origin);
            % Get various sizes needed later
            ncell = model.G.cells.num;
            ncomp = model.getNumberOfComponents();
        end

        function dwdp = getWeightDerivatives(prop, w, w0, factor)
            % Compute weight derivatives by numerical differentiation
            if ~isempty(w0)
                dw = w - w0;
                if prop.fullDerivatives
                    masses = factor.mass;
                    m0 = value(factor.mass0');
                    m  = value(masses');
                    dm = m - m0;
                    dm(abs(dm) < 1e-8) = 0;
                    dwdm = bsxfun(@rdivide, dw, dm);
                    for i = 1:size(dwdm, 2)
                        dwdm(:, i) = dwdm(:, i).*getDiagonal(masses{i}, i);
                    end
                    dwdp = dwdm;
                else
                    p0 = factor.pressure0;
                    p  = factor.pressure;
                    dp = p - p0;
                    dwdp = bsxfun(@rdivide, dw, dp);
                end
                dwdp(~isfinite(dwdp)) = 0;
            else
                dwdp = [];
            end
        end
    end
end

function d = getDiagonal(x, i)
    if issparse(x.jac{1})
        d = full(diag(x.jac{i}));
    else
        d = x.jac{1}.diagonal(:, i);
    end
end

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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

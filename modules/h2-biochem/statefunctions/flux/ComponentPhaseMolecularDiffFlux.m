classdef ComponentPhaseMolecularDiffFlux < StateFunction
    % Molecular diffusion flux for each component/phase.
    %
    % J_{i,alpha} = - K_{i,alpha} * Grad(z_{i,alpha})
    % K_{i,alpha} = (phi*S_alpha*tau_MQ) * rho_alpha * D_{i,alpha}
    %
    % Millington-Quirk:
    %   (phi*S)*tau_MQ = (phi*S)^(7/3) / phi^2

    methods
        function sf = ComponentPhaseMolecularDiffFlux(model, varargin)
            sf@StateFunction(model, varargin{:});

            % Core dependencies
            sf = sf.dependsOn('Density', 'PVTPropertyFunctions');

            % Porosity in h2-biochem is a FlowDiscretization statefunction (BactPorosity)
            % Do not depend on Porosity from FlowDiscretization. Porosity may be static
            % (numeric rock.poro) or dynamic (function handle rock.poro(p, nbact)).

            if isprop(model, 'rock') && isa(model.rock.poro, 'function_handle')
                sf = sf.dependsOn('nbact', 'state');
            end

            % State deps
            sf = sf.dependsOn('s', 'state');
            sf = sf.dependsOn('x', 'state');
            sf = sf.dependsOn('y', 'state');
            sf = sf.dependsOn('pressure', 'state');
            sf = sf.dependsOn('temperature', 'state');

            sf.label = 'J_{i,\\alpha}';
        end

        function J = evaluateOnDomain(sf, model, state)
            if ~isfield(state, 'x')
                ncomp = model.getNumberOfComponents();
                nph   = model.getNumberOfPhases();
                J = cell(ncomp, nph);
                [J{:}] = deal(0);
                return;
            end
            J = localMolecularDiffusionFluxes(sf, model, state);
        end
    end
end

function J = localMolecularDiffusionFluxes(sf, model, state)
    op    = model.operators;
    ncomp = model.getNumberOfComponents();
    nph   = model.getNumberOfPhases();
    nm    = model.getPhaseNames();

    % --- Porosity (dimensionless) from rock (static or dynamic)

    if isprop(model, 'rock') && isa(model.rock.poro, 'function_handle')
        [p, nbact] = model.getProps(state, 'pressure', 'nbact');
        phi = model.rock.poro(p, nbact);
    else
        phi = model.rock.poro;
    end
    phi_safe = max(phi, 1e-12);

    % --- Density
    rho = sf.getEvaluatedExternals(model, state, 'Density');

    % Mole fractions
    [xc, yc] = localGetMoleFractions(model, state);

    % Diffusion parameters
    [Dliq_ref, paramLJ] = localLoadDiffusionDatabase(model);
    Dij_ref = localBinaryDiffusionReference(model, paramLJ);

    % Gas scaling
    [p, T] = model.getProps(state, 'pressure', 'temperature');
    Tref = 273.15 + 40;
    pref = atm;
    p_safe = max(p, 1e-8*barsa);
    gasScale = (T./Tref).^1.75 .* (pref./p_safe);

    % Initialize
    J = cell(ncomp, nph);
    [J{:}] = deal(0);

    L_ix = model.getLiquidIndex();
    V_ix = model.getVaporIndex();

    for ph = 1:nph
        s = model.getProp(state, ['s', nm(ph)]);

        % (phi*S)*tau_MQ
        phiS = phi .* s;
        geom = (phiS).^(7/3) .* (phi_safe).^(-2);

        if ph == L_ix
            for c = 1:ncomp
                Kc = geom .* rho{ph} .* Dliq_ref(c);
                Kf = localFaceHarmonicAvg(op, Kc);
                J{c, ph} = -Kf .* op.Grad(xc{c});
            end

        elseif ph == V_ix
            yAll = model.getProps(state, 'y');
            for c = 1:ncomp
                invDi = 0;
                for j = 1:ncomp
                    if j == c, continue; end
                    yj  = localGetComponentVector(yAll, j);
                    Dij = Dij_ref(c, j) .* gasScale;
                    invDi = invDi + yj ./ max(Dij, 1e-30);
                end
                Di_mix = 1 ./ max(invDi, 1e-30);

                Kc = geom .* rho{ph} .* Di_mix;
                Kf = localFaceHarmonicAvg(op, Kc);
                J{c, ph} = -Kf .* op.Grad(yc{c});
            end
        end
    end
end

function [xc, yc] = localGetMoleFractions(model, state)
    if isfield(state, 'rs') || isfield(state, 'rv')   
        error('Black-oil diffusion with rs/rv is not supported here.');
    end
    [x, y] = model.getProps(state, 'x', 'y');
    if iscell(x)
        xc = x;
    else
        xc = mat2cell(x, size(x, 1), ones(1, size(x, 2)));
    end
    if iscell(y)
        yc = y;
    else
        yc = mat2cell(y, size(y, 1), ones(1, size(y, 2)));
    end
end

function v = localGetComponentVector(z, j)
    if iscell(z)
        v = z{j};
    else
        v = z(:, j);
    end
end

function Kf = localFaceHarmonicAvg(op, Kc)
    if isfield(op, 'N')
        N  = op.N;
        i1 = N(:, 1); i2 = N(:, 2);
        K1 = Kc(i1);  K2 = Kc(i2);
        Kf = 2 .* K1 .* K2 ./ max(K1 + K2, 1e-30);
    else
        Kf = op.faceAvg(Kc);
    end
end

function [Dliq_ref, paramLJ] = localLoadDiffusionDatabase(model)
    namecp = model.compFluid.names();
    ncomp  = numel(namecp);

    Dliq_ref = 1e-9 * ones(ncomp, 1);
    paramLJ  = repmat([3.5, 150.0], ncomp, 1);

    db.Dliq = struct('h2',6.44e-9,'c1',2.15e-9,'methane',2.15e-9,'co2',2.72e-9, ...
        'h2o',3.29e-9,'water',3.29e-9,'n2',2.86e-9,'c2',1.72e-9, ...
        'ethane',1.72e-9,'c3',1.43e-9,'propane',1.43e-9,'h2s',2.15e-9, ...
        'nc4',1.15e-9,'butane',1.15e-9);

    db.LJ = struct('h2',[2.92,59.7],'c1',[3.758,148.6],'methane',[3.758,148.6], ...
        'co2',[3.996,195.2],'h2o',[2.641,809.1],'water',[2.641,809.1], ...
        'n2',[3.798,71.4],'c2',[4.443,215.7],'ethane',[4.443,215.7], ...
        'c3',[5.118,237.1],'propane',[5.118,237.1],'h2s',[3.60,301.0], ...
        'nc4',[5.206,289.5],'butane',[5.206,289.5]);

    for i = 1:ncomp
        key = lower(namecp{i});
        if isfield(db.Dliq, key), Dliq_ref(i) = db.Dliq.(key); end
        if isfield(db.LJ, key),   paramLJ(i,:) = db.LJ.(key);  end
    end
end

function Dij_ref = localBinaryDiffusionReference(model, paramLJ)
    ncomp = size(paramLJ, 1);
    Molmass = 1e3 .* model.compFluid.molarMass; % kg/mol -> g/mol
    sigma = paramLJ(:, 1); % Å

    Dij_ref = zeros(ncomp);
    for i = 1:ncomp
        for j = 1:ncomp
            if i == j
                Dij_ref(i, j) = inf;
                continue;
            end
            sqrtMij   = sqrt(2 * Molmass(i) * Molmass(j) / (Molmass(i) + Molmass(j)));
            sigma_ij2 = 0.25 * (sigma(i) + sigma(j))^2;
            Dij_ref(i, j) = 1e-4 * 0.001858 / (sqrtMij * sigma_ij2);
        end
    end
    Dij_ref = 0.5*(Dij_ref + Dij_ref.');
end

%{
Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.

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
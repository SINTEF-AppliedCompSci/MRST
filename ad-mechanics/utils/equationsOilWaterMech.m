function [eqs, state] = equationsOilWaterMech(p0, sW0, p, sW, pBH, qWs, qOs,  state, model, dt, mechTerm, ...
                                           drivingForces, varargin)

% Equation for oil water system that also takes input from mechanics.

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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


    % Note that state is given only for output
    opt = struct('iteration', -1, ...
                 'resOnly', false); % just to avoid warning
    opt = merge_options(opt, varargin{:});

    W = drivingForces.W;

    s = model.operators;
    G = model.G;
    f = model.fluid;
    rock = model.rock;

    % Evaluate relative permeability
    sO  = 1 - sW;
    sO0 = 1 - sW0;

    [krW, krO] = model.evaluateRelPerm({sW, sO});

    % Multipliers for properties
    [pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

    % Modifiy relperm by mobility multiplier (if any)
    krW = mobMult.*krW; krO = mobMult.*krO;

    % Compute transmissibility
    T = s.T.*transMult;

    % Gravity contribution
    gdz = model.getGravityGradient();

    % Evaluate water properties
    [vW, bW, mobW, rhoW, pW, upcw] = getFluxAndPropsWater_BO(model, p, sW, krW, T, gdz);
    bW0 = model.fluid.bW(p0);

    % Evaluate oil properties
    [vO, bO, mobO, rhoO, p, upco] = getFluxAndPropsOil_BO(model, p, sO, krO, T, gdz);
    bO0 = getbO_BO(model, p0);

    if model.outputFluxes
        state = model.storeFluxes(state, vW, vO, []);
    end
    if model.extraStateOutput
        state = model.storebfactors(state, bW, bO, []);
        state = model.storeMobilities(state, mobW, mobO, []);
        state = model.storeUpstreamIndices(state, upcw, upco, []);
    end

    % EQUATIONS ---------------------------------------------------------------
    % Upstream weight b factors and multiply by interface fluxes to obtain the
    % fluxes at standard conditions.
    bOvO = s.faceUpstr(upco, bO).*vO;
    bWvW = s.faceUpstr(upcw, bW).*vW;

    % Computation of "effective" porosity which take into account the changes
    % due to mechanics.
    effPorVol = rock.poro.*(G.cells.volumes.*pvMult) + rock.alpha .* ...
              mechTerm.new;
    effPorVol0 = rock.poro.*(G.cells.volumes.*pvMult0) + rock.alpha .* ...
        mechTerm.old;

    % Conservation of mass for water
    water = (1./dt).*(effPorVol.*bW.*sW - effPorVol0.*bW0.*sW0)  + s.Div(bWvW);
    % Conservation of mass for oil
    oil = (1./dt).*(effPorVol.*bO.*sO - effPorVol0.*bO0.*sO0)  + s.Div(bOvO);

    eqs = {water, oil};
    
    % Add in any fluxes / source terms prescribed as boundary conditions.
    [eqs, ~, qRes] = addFluxesFromSourcesAndBC(model, eqs, ...
                                               {pW, p},...
                                               {rhoW,     rhoO},...
                                               {mobW,     mobO}, ...
                                               {bW, bO},  ...
                                               {sW, sO}, ...
                                               drivingForces);
    if model.outputFluxes
        state = model.storeBoundaryFluxes(state, qRes{1}, qRes{2}, [], drivingForces);
    end

    % Set up well equations
    if ~isempty(W)
        wm = model.wellmodel;
        wc    = vertcat(W.cells);
        pw   = p(wc);
        rhos = [model.fluid.rhoWS, model.fluid.rhoOS];
        bw   = {bW(wc), bO(wc)};
        mw   = {mobW(wc), mobO(wc)};
        s = {sW(wc), sO(wc)};

        [cqs, weqs, ctrleqs, wc, state.wellSol]  = wm.computeWellFlux(model, ...
                                                          W, state.wellSol, ...
                                                          pBH, {qWs, qOs}, pw, ...
                                                          rhos, bw, mw, s, ...
                                                          {}, 'nonlinearIteration', ...
                                                          opt.iteration);
        % Store the well equations (relate well bottom hole pressures to
        % influx).
        eqs(3:4) = weqs;
        % Store the control equations (trivial equations ensuring that each
        % well will have values corresponding to the prescribed value)
        eqs{5} = ctrleqs;
        % Add source terms to the equations. Negative sign may be
        % surprising if one is used to source terms on the right hand side,
        % but this is the equations on residual form.
        eqs{1}(wc) = eqs{1}(wc) - cqs{1};
        eqs{2}(wc) = eqs{2}(wc) - cqs{2};
    end
end


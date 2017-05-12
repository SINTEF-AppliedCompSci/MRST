function [eqs, state] = equationsBlackOilBiot(state0, st0, p, sW, x,  rs, rv, ...
                                              st, qWell, bhp, wellVars, state, ...
                                              model, dt, mechTerm, drivingForces, ...
                                              varargin);
    
% Equation for black-oil system that also takes input from mechanics.
    
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

    opt = struct('Verbose',     mrstVerbose,...
                 'resOnly',     false,...
                 'iteration',   -1);

    opt = merge_options(opt, varargin{:});

    % Shorter names for some commonly used parts of the model and forces.
    s = model.operators;
    f = model.fluid;
    rock = model.rock;
    G = model.G;
    W = drivingForces.W;

    [p0, sW0, sG0, rs0, rv0, wellSol0] = model.getProps(state0, 'pressure', ...
                                                                'water', 'gas', ...
                                                                'rs', 'rv', 'wellsol');
    
    [sG, rs, rv, rsSat, rvSat] = calculateHydrocarbonsFromStatusBO(model, st, ...
                                                      1-sW, x, rs, rv, p);
    
    % Evaluate relative permeability
    sO  = 1 - sW  - sG;
    sO0 = 1 - sW0 - sG0;
    [krW, krO, krG] = model.evaluateRelPerm({sW, sO, sG});

    % Multipliers for properties
    [pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

    % Modifiy relperm by mobility multiplier (if any)
    krW = mobMult.*krW; krO = mobMult.*krO; krG = mobMult.*krG;

    % Compute transmissibility
    T = s.T.*transMult;

    % Gravity gradient per face
    gdz = model.getGravityGradient();

    % Evaluate water properties
    [vW, bW, mobW, rhoW, pW, upcw] = getFluxAndPropsWater_BO(model, p, sW, krW, T, gdz);
    bW0 = f.bW(p0);

    % Evaluate oil properties
    [vO, bO, mobO, rhoO, p, upco] = getFluxAndPropsOil_BO(model, p, sO, krO, T, gdz, rs, ~st{1});
    bO0 = getbO_BO(model, p0, rs0, ~st0{1});

    % Evaluate gas properties
    bG0 = getbG_BO(model, p0, rv0, ~st0{2});
    [vG, bG, mobG, rhoG, pG, upcg] = getFluxAndPropsGas_BO(model, p, sG, krG, T, gdz, rv, ~st{2});

    % Store fluxes / properties for debugging / plotting, if requested.
    if model.outputFluxes
        state = model.storeFluxes(state, vW, vO, vG);
    end
    if model.extraStateOutput
        state = model.storebfactors(state, bW, bO, bG);
        state = model.storeMobilities(state, mobW, mobO, mobG);
        state = model.storeUpstreamIndices(state, upcw, upco, upcg);
    end

    % EQUATIONS -----------------------------------------------------------

    % Upstream weight b factors and multiply by interface fluxes to obtain the
    % fluxes at standard conditions.
    bOvO = s.faceUpstr(upco, bO).*vO;
    bWvW = s.faceUpstr(upcw, bW).*vW;
    bGvG = s.faceUpstr(upcg, bG).*vG;


    % Computation of "effective" porosity which take into account the changes
    % due to mechanics.
    effPorVol = rock.poro.*(G.cells.volumes.*pvMult) + rock.alpha .* ...
              mechTerm.new;
    effPorVol0 = rock.poro.*(G.cells.volumes.*pvMult0) + rock.alpha .* ...
        mechTerm.old;
    
    % The first equation is the conservation of the water phase. This equation is
    % straightforward, as water is assumed to remain in the aqua phase in the
    % black oil model.

    water = (1./dt).*(effPorVol.*bW.*sW - effPorVol0.*bW0.*sW0)  + s.Div(bWvW);

    % Second equation: mass conservation equation for the oil phase at surface
    % conditions. This is any liquid oil at reservoir conditions, as well as
    % any oil dissolved into the gas phase (if the model has vapoil enabled).
    if model.vapoil
        % The model allows oil to vaporize into the gas phase. The conservation
        % equation for oil must then include the fraction present in the gas
        % phase.
        rvbGvG = s.faceUpstr(upcg, rv).*bGvG;
        % Final equation
        oil = (1./dt).*( effPorVol.* (bO.* sO  + rv.* bG.* sG) - ...
                           effPorVol0.*(bO0.*sO0 + rv0.*bG0.*sG0) ) + ...
              s.Div(bOvO + rvbGvG);
    else
        oil = (1./dt).*( effPorVol.*bO.*sO - effPorVol0.*bO0.*sO0 ) + s.Div(bOvO);
    end

    % Conservation of mass for gas. Again, we have two cases depending on
    % whether the model allows us to dissolve the gas phase into the oil phase.
    if model.disgas
        % The gas transported in the oil phase.
        rsbOvO = s.faceUpstr(upco, rs).*bOvO;
        
        gas = (1./dt).*( effPorVol.* (bG.* sG  + rs.* bO.* sO) - ...
                           effPorVol0.*(bG0.*sG0 + rs0.*bO0.*sO0 ) ) + ...
              s.Div(bGvG + rsbOvO);
    else
        gas = (1./dt).*( effPorVol.*bG.*sG - effPorVol0.*bG0.*sG0 ) + s.Div(bGvG);
    end

    % Put the set of equations into cell arrays along with their names/types.
    eqs = {water, oil, gas};
    names = {'water', 'oil', 'gas'};
    types = {'cell', 'cell', 'cell'};

    % Add in any fluxes / source terms prescribed as boundary conditions.
    [eqs, ~, qRes] = addFluxesFromSourcesAndBC(model, eqs, ...
                                               {pW, p, pG},...
                                               {rhoW,     rhoO, rhoG},...
                                               {mobW,     mobO, mobG}, ...
                                               {sW, sO, sG}, ...
                                               drivingForces);
    if model.outputFluxes
        state = model.storeBoundaryFluxes(state, qRes{1}, qRes{2}, qRes{3}, drivingForces);
    end

    % Finally, add in and setup well equations
    if ~isempty(W)
        dissolved = model.getDissolutionMatrix(rs, rv);
        wellSol = model.getProp(state, 'wellsol');
        [~, ~, wellVars, wellVarNames, wellMap] = ...
            model.FacilityModel.getAllPrimaryVariables(wellSol);
        [eqs, ~, ~, state.wellSol] = model.insertWellEquations(eqs, names, ...
                                                          types, wellSol0, ...
                                                          wellSol, qWell, bhp, ...
                                                          wellVars, wellMap, ...
                                                          p, {mobW, mobO, ...
                            mobG}, {rhoW, rhoO, rhoG}, dissolved, {}, dt, opt);
    end
    
end


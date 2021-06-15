function [eqs, names, types, state] = equationsOilWaterMech(p0, sW0, state0, ...
                                                            p, sW, wellVars, state, ...
                                                            model, dt, mechTerm, ...
                                                            drivingForces, varargin)
%
%
% SYNOPSIS:
%   function [eqs, names, types, state] = equationsOilWaterMech(p0, sW0, state0, p, sW, wellVars, state, model, dt, mechTerm, drivingForces, varargin)
%
% PARAMETERS:
%   p0            - pressure   (for previous time step)  
%   sW0           - saturation (for previous time step)
%   state0        - state      (for previous time step)
%   p             - pressure
%   sW            - saturation
%   wellVars      - well variables
%   state         - current state
%   model         - model class instance that is used
%   dt            - time step size
%   mechTerm      - mechanical input which will enter the computation of the
%                   effective porosity
%   drivingForces - structure that gathers the well parameters and boundary conditions.
%   varargin      - 
%
% RETURNS:
%   eqs   - The residual values as ADI variables (that is with the Jacobian)
%           if the inputs were also ADI.
%   names - The name of each equations
%   types - The type of each equations
%   state - Some field related to well control of the state variables may be updated.
%
% EXAMPLE:
%   run2DCase, runNorneExample
%
% SEE ALSO:
%   equationOilWater.

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
    names = {'water', 'oil'};
    types = {'cell', 'cell'};
    
    % Finally, add in and setup well equations
    wellSol = model.getProp(state, 'wellsol');
    [~, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);
    wellSol0 = model.getProp(state0, 'wellsol');
    
    [eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, ...
                                                      types, wellSol0, wellSol, ...
                                                      wellVars, wellMap, p, ...
                                                      {mobW, mobO}, {rhoW, ...
                        rhoO}, {}, {}, dt, opt);
    
    
end


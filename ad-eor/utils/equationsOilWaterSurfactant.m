function [problem, state] = equationsOilWaterSurfactant(state0, state, model, ...
                                                      dt, drivingForces, ...
                                                      varargin)
%
% SYNOPSIS:
%   function [problem, state] = equationsOilWaterSurfactant(state0, state, model, dt, drivingForces, varargin)
%
% DESCRIPTION: Assemble the linearized equations for an oil-water-surfactant
% system, computing both the residuals and the Jacobians. Returns the result as
% an instance of the class LinearizedProblem which can be solved using instances
% of LinearSolverAD.
%
% A description of the modeling equations can be found in the directory
% ad-eor/docs.
%
%
% PARAMETERS:
%   state0        - State at previous times-step
%   state         - State at current time-step
%   model         - Model instance
%   dt            - time-step
%   drivingForces - Driving forces (boundary conditions, wells, ...)
%   varargin      - optional parameters
%
% RETURNS:
%   problem - Instance of LinearizedProblem
%   state   - Updated state variable (fluxes, mobilities and more can be
%             stored, the wellSol structure is also updated in case of control switching)
%
% EXAMPLE:
%
% SEE ALSO: LinearizedProblem, LinearSolverAD
%
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

    opt = struct('Verbose', mrstVerbose, ...
                 'reverseMode', false, ...
                 'resOnly', false, ...
                 'iteration', -1 );
    opt = merge_options(opt, varargin{:});

    W     = drivingForces.W;
    fluid = model.fluid;
    op    = model.operators;
    G     = model.G;

    % Properties at current timestep
    [p, sW, c, cmax, wellSol] = model.getProps(state, 'pressure', 'water', ...
                                                      'surfactant', ...
                                                      'surfactantmax', ...
                                                      'wellsol');

    % Properties at previous timestep
    [p0, sW0, c0, cmax0, wellSol0] = model.getProps(state0, 'pressure', 'water', ...
                                                            'surfactant', ...
                                                            'surfactantmax', ...
                                                            'wellsol');

    [qWell, pBH, wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);
    % Typically the primary well variables are :
    %  - the phase well rates (qWell)
    %  - the bottom hole pressures
    %  - the surfactant concentrations, at injection and production wells,
    %    contained in the wellVars, wellVarNames, wellMap structures

    % Initialize independent variables.
    if ~opt.resOnly,
        % ADI variables needed since we are not only computing residuals.
        if ~opt.reverseMode,
            [p, sW, c, qWell{:}, pBH, wellVars{:}] = initVariablesADI(p, sW, ...
                                                              c,  qWell{:}, ...
                                                              pBH, wellVars{:});
        else
            % Not implemented yet
            % zw = zeros(size(pBH));
            % [p0, sW0, c0, zw, zw, zw, zw] = ...
            %     initVariablesADI(p0, sW0, c0, zw, zw, zw, zw); %#ok
            % clear zw
        end
    end

    % We will solve for pressure, water saturation (oil saturation follows via
    % the definition of saturations), surfactant concentration and well rates +
    % bhp.
    primaryVars = {'pressure', 'sW', 'surfactant', wellVarNames{:}};



    % EQUATIONS ---------------------------------------------------------------

    % Compute fluxes and other properties for oil and water.
    [dp, mob, upc, b, rho, pvMult, b0, pvMult0, T] = ...
        computeFluxAndPropsOilWaterSurfactant(model, p0, p, sW, c, pBH, W);

    dpW  = dp{1} ; dpO  = dp{2};
    mobW = mob{1}; mobO = mob{2};
    rhoW = rho{1}; rhoO = rho{2};
    upcW = upc{1}; upcO = upc{2};
    bW   = b{1}  ; bO   = b{2};
    bW0  = b0{1} ; bO0  = b0{2};


    % Upstream weight b factors and multiply by interface fluxes to obtain the
    % fluxes at standard conditions.
    vO     =  -op.faceUpstr(upcO, mobO).*T.*dpO;
    vW     = -op.faceUpstr(upcW, mobW).*T.*dpW;
    bOvO   = op.faceUpstr(upcO, bO).*vO;
    bWvW   = op.faceUpstr(upcW, bW).*vW;

    % Conservation of mass for water
    water = (op.pv/dt).*(pvMult.*bW.*sW - pvMult0.*bW0.*sW0) + op.Div(bWvW);

    % Conservation of mass for oil
    sO  = 1 - sW;
    sO0 = 1 - sW0;
    oil = (op.pv/dt).*(pvMult.*bO.*sO - pvMult0.*bO0.*sO0) + op.Div(bOvO);

    % Computation of adsoprtion term
    poro = model.rock.poro;
    ads  = effads(c, cmax, fluid);
    ads0 = effads(c0, cmax0, fluid);
    ads_term = fluid.rhoRSft.*((1-poro)./poro).*(ads - ads0);

    % Conservation of surfactant in water:
    mobSft = mobW.*c;
    vSft   = -op.faceUpstr(upcW, mobSft).*T.*dpW;
    bWvSft = op.faceUpstr(upcW, bW).*vSft;
    surfactant = (op.pv/dt).*((pvMult.*bW.*sW.*c - pvMult0.*bW0.*sW0.*c0) +  ads_term) + ...
        op.Div(bWvSft);

    if model.extraStateOutput
        sigma = fluid.ift(c);
    end

    eqs   = {water, oil, surfactant};
    names = {'water', 'oil', 'surfactant'};
    types = {'cell', 'cell', 'cell'};

    % Finally, add in and setup well equations
    if ~isempty(W)
        [eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, ...
                                                          names, types, wellSol0, ...
                                                          wellSol, qWell, pBH, ...
                                                          wellVars, wellMap, ...
                                                          p, mob, rho, {}, ...
                                                          {c}, dt, opt);
    end
    problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

end

% Effective adsorption, depending of desorption or not
function y = effads(c, cmax, fluid)
   if fluid.adsInxSft == 2
      y = fluid.surfads(max(c, cmax));
   else
      y = fluid.surfads(c);
   end
end

function [problem, state] = equationsOilWaterSurfactant(state0, state, model, ...
                                                      dt, drivingForces, ...
                                                      varargin)
%
% SYNOPSIS:
%   function [problem, state] = equationsOilWaterSurfactant(state0, state, model, dt, drivingForces, varargin)
%
% DESCRIPTION: 
%   Assemble the linearized equations for an oil-water-surfactant
%   system, computing both the residuals and the Jacobians. Returns the result as
%   an instance of the class LinearizedProblem which can be solved using instances
%   of LinearSolverAD.
%
%   A description of the modeling equations can be found in the directory
%   ad-eor/docs.
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

    opt = struct('Verbose'        , mrstVerbose, ...
                 'reverseMode'    , false      , ...
                 'velocCompMethod', 'square'   , ...
                 'resOnly'        , false      , ...
                 'iteration'      , -1 );
    opt = merge_options(opt, varargin{:});

    G     = model.G;
    op    = model.operators;
    fluid = model.fluid;
    W     = drivingForces.W;

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

    [wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);
    % Typically the primary well variables are :
    %  - the phase well rates (qWell)
    %  - the bottom hole pressures
    %  - the surfactant concentrations, at injection and production wells,
    %    contained in the wellVars, wellVarNames, wellMap structures

    % Initialize independent variables.
    if ~opt.resOnly,
        % ADI variables needed since we are not only computing residuals.
        if ~opt.reverseMode,
            [p, sW, c, wellVars{:}] = initVariablesADI(p, sW, c, wellVars{:});
        else
            wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
            [p0, sW0, c0, wellVars0{:}] = initVariablesADI(p0, sW0, c0, wellVars0{:}); %#ok
        end
    end

    % We will solve for pressure, water saturation (oil saturation follows via
    % the definition of saturations), surfactant concentration and well rates +
    % bhp.
    primaryVars = {'pressure', 'sW', 'surfactant', wellVarNames{:}};

    sO  = 1 - sW;
    sO0 = 1 - sW0;
    sat  = {sW, sO};
    sat0 = {sW0, sO0};
    
    % Update state with AD-variables
    state = model.setProps(state  , {'s', 'pressure', 'surfactant'}, {sat , p , c});
    state0 = model.setProps(state0, {'s', 'pressure', 'surfactant'}, {sat0, p0, c0});
    % Set up properties
    state = model.initStateFunctionContainers(state);
    
    % EQUATIONS ---------------------------------------------------------------
    pBH = wellVars{wellMap.isBHP};
    Nc = computeCapillaryNumber(p, c, pBH, W, fluid, G, op, 'velocCompMethod', ...
                                opt.velocCompMethod);
    state.CapillaryNumber = Nc;

    [b, pv] = model.getProps(state, 'ShrinkageFactors', 'PoreVolume');
    [b0, pv0] = model.getProps(state0, 'ShrinkageFactors', 'PoreVolume');
    [phaseFlux, flags] = model.getProps(state, 'PhaseFlux',  'PhaseUpwindFlag');
    [pressures, mob, rho] = model.getProps(state, 'PhasePressures', 'Mobility', 'Density');

    [bW, bO]     = deal(b{:});
    [bW0, bO0]   = deal(b0{:});
    [vW, vO]     = deal(phaseFlux{:});
    [upcw, upco] = deal(flags{:});
    [mobW, mobO] = deal(mob{:});
    
    % Upstream weight b factors and multiply by interface fluxes to obtain the
    % fluxes at standard conditions.
    bOvO   = op.faceUpstr(upco, bO).*vO;
    bWvW   = op.faceUpstr(upcw, bW).*vW;

    % Conservation of mass for water
    water = (1/dt).*( pv.*bW.*sW - pv0.*bW0.*sW0 );
    divWater = op.Div(bWvW);
    
    % Conservation of mass for oil
    oil = (1/dt).*( pv.*bO.*sO - pv0.*bO0.*sO0 );
    divOil = op.Div(bOvO);
    
    % Computation of adsoprtion term
    poro = model.rock.poro;
    ads  = model.getProp(state , 'SurfactantAdsorption');
    ads0 = model.getProp(state0, 'SurfactantAdsorption');
    ads_term = fluid.rhoRSft.*((1-poro)./poro).*(ads - ads0);

    % Conservation of surfactant in water:
    vSft   = op.faceUpstr(upcw, c).*vW;
    bWvSft = op.faceUpstr(upcw, bW).*vSft;
    surfactant = (1/dt)*((pv.*bW.*sW.*c - pv0.*bW0.*sW0.*c0) + ads_term);
    divSurfactant = op.Div(bWvSft);
    
    if model.extraStateOutput
        sigma = fluid.ift(c);
    end

    eqs        = {water   , oil   , surfactant};
    divTerms   = {divWater, divOil, divSurfactant};
    names      = {'water' , 'oil' , 'surfactant'};
    types      = {'cell'  , 'cell', 'cell'};
    components = {c};
    
    [eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, ...
                                                   state, pressures, sat, mob, ...
                                                   rho, {}, components, ...
                                                   drivingForces);
    % Finally, add in and setup well equations
    [eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, ...
                                                      types, wellSol0, wellSol, ...
                                                      wellVars, wellMap, p, ...
                                                      mob, rho, {}, components, dt, ...
                                                      opt);
    
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

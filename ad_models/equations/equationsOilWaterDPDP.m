function [problem, state] = equationsOilWaterDPDP(state0, state, model, dt, drivingForces, varargin)
% Generate linearized problem for the two-phase oil-water dual-porosity dual permeability model
%
% SYNOPSIS:
%   [problem, state] = equationsOilWater(state0, state, model, dt, drivingForces)
%
% DESCRIPTION:
%   This is the core function of the two-phase oil-water solver. This
%   function assembles the residual equations for the conservation of water
%   and oil as well as required well equations. By default, Jacobians are
%   also provided by the use of automatic differentiation.
%
% REQUIRED PARAMETERS:
%   state0    - Reservoir state at the previous timestep. Assumed to have
%               physically reasonable values.
%
%   state     - State at the current nonlinear iteration. The values do not
%               need to be physically reasonable.
%
%   model     - TwoPhaseOilWaterModel-derived class. Typically,
%               equationsOilWater will be called from the class
%               getEquations member function.
%
%   dt        - Scalar timestep in seconds.
%
%   drivingForces - Struct with fields:
%                   * W for wells. Can be empty for no wells.
%                   * bc for boundary conditions. Can be empty for no bc.
%                   * src for source terms. Can be empty for no sources.
%
% OPTIONAL PARAMETERS:
%   'Verbose'    -  Extra output if requested.
%
%   'reverseMode'- Boolean indicating if we are in reverse mode, i.e.
%                  solving the adjoint equations. Defaults to false.
%
%   'resOnly'    - Only assemble residual equations, do not assemble the
%                  Jacobians. Can save some assembly time if only the
%                  values are required.
%
%   'iterations' - Nonlinear iteration number. Special logic happens in the
%                  wells if it is the first iteration.
% RETURNS:
%   problem - LinearizedProblemAD class instance, containing the water
%             and oil conservation equations, as well as well equations
%             specified by the WellModel class.
%
%   state   - Updated state. Primarily returned to handle changing well
%             controls from the well model.
%
% SEE ALSO:
%   equationsBlackOil, TwoPhaseOilWaterModel


opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'iteration', -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.W;
s = model.operators;
mo = model.operators_matrix;
G = model.G;
f = model.fluid_matrix; % f is used for obtaining matrix properties


% Properties at current timestep
[p, sW, wellSol] = model.getProps(state, 'pressure', 'water', 'wellsol');
% Properties at previous timestep
[p0, sW0, wellSol0] = model.getProps(state0, 'pressure', 'water', 'wellsol');

% Set up current time in the states
if ~isfield(state, 'time')
    state0.time = 0;
    state.time = dt;
else
    if isfield(state0, 'time')
        state.time = state0.time + dt;
    end
    % else we are within the Newton's iterations
end

% Matrix properties
[pom,swm] = model.getProps(state, 'pm','swm');
[pom0,swm0] = model.getProps(state0, 'pm','swm');

[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

% Initialize independent variables.
if true  %~opt.resOnly
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode
        [p, sW, pom, swm, wellVars{:}] = initVariablesADI(p, sW, pom, swm, wellVars{:});
    else
        wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
        [p0, sW0, pom0, swm0, wellVars0{:}] = initVariablesADI(p0, sW0, pom0, swm0, wellVars0{:}); %#ok
        clear zw
    end
else
    msg = 'Maximum number of iterations reached... Exiting..';
    %error(msg)
end

% We will solve for pressure, water saturation (oil saturation follows via
% the definition of saturations) and well rates + bhp.
primaryVars = {'pressure', 'sW', 'pm', 'swm', wellVarNames{:}};

%% Properties for Fracture
% Evaluate relative permeability
sO  = 1 - sW;
sO0 = 1 - sW0;

[krW, krO] = model.evaluateRelPerm({sW, sO},'medium','fracture');

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

%% Properties for Matrix
% Evaluate relative permeability
som  = 1 - swm;
som0 = 1 - swm0;

[krWm, krOm] = model.evaluateRelPerm({swm, som},'medium','matrix');

% Multipliers for properties
[pvMultm, transMultm, mobMultm, pvMultm0] = getMultipliers(model.fluid_matrix, pom, pom0);

% Modifiy relperm by mobility multiplier (if any)
krWm = mobMultm.*krWm; krOm = mobMultm.*krOm;

% Compute transmissibility
Tm = mo.T.*transMultm;

% Gravity contribution
gdzm = model.getGravityGradient();

% Evaluate water properties
[vWm, bWm, mobWm, rhoWm, pwm, upcwm] = getFluxAndPropsWater_BO(model, pom, swm, krWm, Tm, gdzm);
bWm0 = model.fluid_matrix.bW(pom0);

% Evaluate oil properties
[vOm, bOm, mobOm, rhoOm, pom, upcom] = getFluxAndPropsOil_BO(model, pom, som, krOm, Tm, gdzm);
bOm0 = getbO_BO(model, pom0);

%% Transfer
vb = model.G.cells.volumes;

matrix_fields.pom = pom;
matrix_fields.swm = swm;
fracture_fields.pof = p;
fracture_fields.swf = sW;

transfer_model = model.transfer_model_object;

[Talpha] = transfer_model.calculate_transfer(model,fracture_fields,matrix_fields);

Twm = vb.*Talpha{1};
Tom = vb.*Talpha{2};

%% Output Additional Info
if model.outputFluxes
    
    if ~isempty(drivingForces.bc)
    % Since the parent class ReservoirModel only knows of state.flux
    % (interpreted as the fracture flux here), we explicitly save matrix 
    % fluxes in state.flux_matrix
    
    % Save internal fluxes
    state = model.storeFluxes(state, vWm, vOm, []);
    state.flux_matrix = state.flux;
    state = model.storeFluxes(state, vW, vO, []);
    
    % Save boundary fluxes from pressure gradients
    % Temporarily set the operators for the matrix
    frop = model.operators;
    model.operators = model.operators_matrix;
    [~, ~, ~, qRes] = getBoundaryConditionFluxesAD(model, ...
        {pwm, pom}, {swm, som}, {mobWm, mobOm}, {rhoWm, rhoOm}, {bWm, bOm}, ...
        drivingForces.bc);
    model.operators = frop; % Restore the fracture operators
    state = model.storeBoundaryFluxes(state, qRes{1}, qRes{2}, [], drivingForces);
    state.flux_matrix = state.flux;
    [~, ~, ~, qRes] = getBoundaryConditionFluxesAD(model, ...
        {pW, p}, {sW, sO}, {mobW, mobO}, {rhoW, rhoO}, {bW, bO}, ...
        drivingForces.bc);
    state = model.storeBoundaryFluxes(state, qRes{1}, qRes{2}, [], drivingForces);
    
    end

end

if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, []);
    state = model.storeMobilities(state, mobW, mobO, []);
    state = model.storeUpstreamIndices(state, upcw, upco, []);
    state.Twm = value(Twm);
    state.Tom = value(Tom);
    %state.pcm = value(pcOWm); % Do not consider pc for the moment
    state.pW = value(pW);
    state.pwm = value(pwm);
end

%% EQUATIONS ---------------------------------------------------------------
% Fracture
% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.

% Fracture
bOvO = s.faceUpstr(upco, bO).*vO;
bWvW = s.faceUpstr(upcw, bW).*vW;
% Matrix
bOvOm = s.faceUpstr(upcom, bOm).*vOm;
bWvWm = s.faceUpstr(upcwm, bWm).*vWm;

% Conservation of mass for water - fracture
water_fracture = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + s.Div(bWvW);
water_fracture = water_fracture + Twm;

% Conservation of mass for oil - fracture
oil_fracture = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) + s.Div(bOvO);
oil_fracture = oil_fracture + Tom;

eqs{1} = water_fracture;
eqs{2} = oil_fracture;

% Matrix
% Conservation of mass for water - matrix
water_matrix = (s.pv_matrix/dt).*( pvMultm.*bWm.*swm - pvMultm0.*bWm0.*swm0 ) + mo.Div(bWvWm);
water_matrix = water_matrix - Twm;

% Conservation of mass for oil - matrix
oil_matrix = (s.pv_matrix/dt).*( pvMultm.*bOm.*som - pvMultm0.*bOm0.*som0 ) + mo.Div(bOvOm);
oil_matrix = oil_matrix - Tom;

eqs{3} = water_matrix;
eqs{4} = oil_matrix;

names = {'water', 'oil','water_matrix', 'oil_matrix'};
types = {'cell', 'cell', 'cell', 'cell'};

% Add in any fluxes / source terms prescribed as boundary conditions.
rho = {rhoW, rhoO};
mob = {mobW, mobO};
sat = {sW, sO};
          

% Compute the contribution from BC whil avoiding the capillary jump
pcOW = 0;
if isfield(model.fluid, 'pcOW')
    pcOW  = model.fluid.pcOW(sW);
end

% Include BC and sources in fracture equations 
[eqsfrac, state] = addBoundaryConditionsAndSources(model, ...
                eqs(1:2), names(1:2), types(1:2), state, ...                                                                 
                {pW + pcOW, p}, {sW, sO}, {mobW, mobO}, {rhoW, rhoO}, ...    
                {}, {}, ...
                drivingForces);                                                           

eqs(1:2) = eqsfrac;

if ~isempty(drivingForces.bc)

    % Compute the contribution from BC while avoiding the capillary jump
    pcOW = 0;
    if isfield(model.fluid, 'pcOW')
        pcOW  = model.fluid.pcOW(swm);
    end

    % Store the fracture operators and BC values
    frop = model.operators;
    frbc = drivingForces.bc.value;

    % Set up the boundary conditions using the matrix operators and BC values
    model.operators = model.operators_matrix;
    drivingForces.bc.value = drivingForces.bc.value_matrix;

    [eqsmat, state] = addBoundaryConditionsAndSources(model, ...
                    eqs(3:4), {'water', 'oil'}, types(3:4), state, ...                                                              
                    {pwm + pcOW, pom}, {swm, som}, {mobWm, mobOm}, {rhoWm, rhoOm}, ...    
                    {}, {}, ...
                    drivingForces);                                                           

    eqs(3:4) = eqsmat;

    % Restore the fracture operators and BC values
    model.operators = frop;
    drivingForces.bc.value = frbc;

end

% Finally, add in and setup well equations
[eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, types, wellSol0, wellSol, wellVars, wellMap, p, mob, rho, {}, {}, dt, opt);

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

end

%{
Copyright 2022 Geological Survey of Denmark and Greenland (GEUS).

Author: Nikolai Andrianov, nia@geus.dk.

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
 
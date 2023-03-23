function [eqs, names, types, state] = equationsDCWaterMech(state0, state, model, ...
                                                           dt, drivingForces, varargin)
%
%
% SYNOPSIS:
%   function [eqs, names, types, state] = equationsDPWaterMech(state0, state, model, dt, drivingForces, varargin)
%
% DESCRIPTION: 
%   Assembles the residuals for the poromechanically coupled dual-permeability mass balance.
%
% PARAMETERS:
%   state0                   - state    (previous time step)
%   state                    - state (current time step)
%   model                    - model class instance that is used.
%                              fluidModel in this case.
%   dt                       - time step
%   drivingForces            - structure that gathers the well parameters and boundary conditions.
%   varargin                 - 
%
% RETURNS:
%   eqs   - The residual values as ADI variables (that is with the Jacobian)
%           if the inputs were also ADI.
%   names - The name of each equations
%   types - The type of each equations
%   state - Some field related to well control of the state variables may be updated.
%
% EXAMPLE: example_void_fractures
%
% SEE ALSO: 
%
%{
Copyright 2009-2020 SINTEF ICT, Applied Mathematics.

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
opt = struct('Verbose', mrstVerbose, ...
         'reverseMode', false,...
         'resOnly', false,...
         'iteration', -1);  % Compatibility only
opt = merge_options(opt, varargin{:});

% Shorter names for some commonly used parts of the model and forces.
f = model.fluidModel.fluid;
s = model.fluidModel.operators;
f_m = model.fluidModel.fluid_matrix;
s_m = model.fluidModel.operators_matrix;    
transferModel = model.transfer_model_object;


%% -------------------------------------------------------------------------
% Initialization of independent variables
% properties at current timestep, and current iteration level
[p, pm, wellSol, xd] = model.getProps(state, 'pressure', 'pressure_matrix', 'wellSol', 'xd');
% properties at previous timestep
[p0, pm0, wellSol0, xd0] = model.getProps(state0, 'pressure', 'pressure_matrix', 'wellSol', 'xd');

% well props
%[frac_index, mat_index] = model.fluidModel.findDCWells(wellSol);
%[wellVars, ~, wellMap] = model.fluidModel.FacilityModel.getAllPrimaryVariables(wellSol(:,frac_index));
%mat_index = 2;
[wellVars, ~, wellMap] = ...
    model.fluidModel.FacilityModel.getAllPrimaryVariables(wellSol);
if ~opt.reverseMode
    [p, pm, wellVars{:}, xd] = initVariablesADI(p, pm, wellVars{:}, xd);
else
    error('Reverse mode AD currently not implemented for DC-mech module')
end

gdz = s.Grad(model.G.cells.centroids) * model.fluidModel.getGravityVector()';

% update state with AD-variables (required for state functions)
state = model.setProps(state, {'pressure', 'pressure_matrix', 'xd'}, {p, pm, xd});
state0 = model.setProps(state0, {'pressure', 'pressure_matrix', 'xd'}, {p0, pm0, xd0});
% set up properties
state = model.initStateFunctionContainers(state);


%% -------------------------------------------------------------------------
% Rock props
[poro_m, poro_f] = model.getProps(state, 'PoroelasticMatrixPoro', 'PoroelasticFracturePoro');
[poro_m0, poro_f0] = model.getProps(state0, 'PoroelasticMatrixPoro', 'PoroelasticFracturePoro');

% check for fracture multipliers:
trMult = 1;
if isfield(f, 'tranMultR'), trMult = f.tranMultR(p); end

transMult=1;
if isfield(f, 'transMult')
   transMult = f.transMult(p);
elseif isfield(model.rock, 'nonLinearPerm')
   transMult = model.getProps(state, 'PoroelasticFractureTransMult');
end
T_f = s.T.*transMult;

% check for matrix multipliers:
trMult_m = 1;
if isfield(f_m, 'tranMultR'), trMult_m = f_m.tranMultR(pm); end


%% -------------------------------------------------------------------------
% Fluid props
% evaluate fracture-water properties
bW     = f.bW(p);
bW0 = f.bW(p0);
rhoW   = bW.*f.rhoWS;
% rhoW on face, avarge of neighboring cells (E100, not E300)
rhoWf  = s.faceAvg(rhoW);
mobW   = trMult./f.muW(p);
dpW     = s.Grad(p) - rhoWf.*gdz;
% water upstream-index
upcw = (value(dpW)<=0);
vW = - s.faceUpstr(upcw, mobW).*T_f.*dpW;
bWvW = s.faceUpstr(upcw, bW).*vW;   

% evaluate matrix-water properties
bWm     = f_m.bW(pm);
bWm0 = f_m.bW(pm0);
rhoWm   = bWm.*f_m.rhoWS;
% rhoW on face, avarge of neighboring cells (E100, not E300)
rhoWfm  = s_m.faceAvg(rhoWm);
mobWm   = trMult_m./f_m.muW(pm);
dpWm     = s_m.Grad(pm) - rhoWfm.*gdz;
% water upstream-index
upcwm = (value(dpWm)<=0);
vWm = - s_m.faceUpstr(upcwm, mobWm).*s_m.T.*dpWm;
bWvWm = s_m.faceUpstr(upcwm, bWm).*vWm;    


%% -------------------------------------------------------------------------
% Mass Transfer model
matrix_fields.pm = pm;
fracture_fields.pf = p;
[Talpha] = transferModel.calculate_transfer(model, fracture_fields, matrix_fields);
Talpha = model.G.cells.volumes.*Talpha;

% Output Additional Information
if model.outputFluxes
    state = model.storeFluxes(state, vW, vWm, []);
end

if model.extraStateOutput
    state = model.storebfactors(state, bW, bWm, []);
    state = model.storeMobilities(state, mobW, mobWm, []);
    state = model.storeUpstreamIndices(state, upcw, upcwm, []);
    state.Talpha = Talpha;
end


%% EQUATIONS ---------------------------------------------------------------
names = {'water', 'water_matrix'};
types = {'cell', 'cell'};
mob = {mobW, mobWm}; 
rho = {rhoW, rhoWm}; 

% Fracture
% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.
eqs{1} = (model.G.cells.volumes./dt).*( poro_f.*bW - poro_f0.*bW0 ) + s.Div(bWvW) ...
          + Talpha;    

% Matrix
% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.
eqs{2} = (model.G.cells.volumes./dt).*( poro_m.*bWm - poro_m0.*bWm0 ) + s.Div(bWvWm) ...
          - Talpha;

% Add BCs
% dummy saturation
sW = ones(model.G.cells.num, 1);
[eqs(1), state] = addBoundaryConditionsAndSources(model.fluidModel, ...
                                                  eqs(1), names(1), types(1), state, ...
                                                  {p}, {sW}, mob(1), rho(1), ...
                                                  {}, {}, ...
                                                  drivingForces);
[eqs(2), state] = addBoundaryConditionsAndSources(model.fluidModel, ...
                                                  eqs(2), names(2), types(2), state, ...
                                                  {pm}, {sW}, mob(2), rho(2), ...
                                                  {}, {}, ...
                                                  drivingForces);

% Well equations, for now, just implemented from fracture  
% these need to come from the fluid model!
[eqs, names, types, state.wellSol] = model.fluidModel.insertWellEquations(...
                                                          eqs, names, types, wellSol0,...
                                                          wellSol, wellVars, wellMap,...
                                                          p, mob(1), rho(1), {}, {}, dt, opt);


end

function [problem, state] = equationsWaterDP(state0, state, model, dt, drivingForces, varargin)
% Generate linearized problem for the single-phase water model
%
% SYNOPSIS:
%   [problem, state] = equationsWater(state0, state, model, dt, drivingForces)
%
% DESCRIPTION:
%   This is the core function of the single-phase water solver with
%   black-oil style properties. This function assembles the residual
%   equations for the conservation of water and oil as well as required
%   well equations. By default, Jacobians are also provided by the use of
%   automatic differentiation.
%
% REQUIRED PARAMETERS:
%   state0    - Reservoir state at the previous timestep. Assumed to have
%               physically reasonable values.
%
%   state     - State at the current nonlinear iteration. The values do not
%               need to be physically reasonable.
%
%   model     - WaterModel-derived class. Typically,
%               equationsWater will be called from the class
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
%   problem - LinearizedProblemAD class instance, containing the equation
%               for the water pressure, as well as well equations specified
%               by the WellModel class.
%
%   state   - Updated state. Primarily returned to handle changing well
%             controls from the well model.
%
% SEE ALSO:
%   equationsBlackOil, ThreePhaseBlackOilModel

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

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'iteration', -1);  % Compatibility only

opt = merge_options(opt, varargin{:});

W = drivingForces.W;

s = model.operators;
G = model.G;
f = model.fluid;
fm = model.fluid_matrix;

%% Dual porosity: adding matrix
[p, wellSol, pm] = model.getProps(state, 'pressure', 'wellsol', 'pressure_matrix');
[p0, wellSol0, pm0] = model.getProps(state0, 'pressure', 'wellSol', 'pressure_matrix');

%%
[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

%Initialization of independent variables ----------------------------------

if ~opt.resOnly
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode
        % Dual porosity: adding matrix
        [p, pm, wellVars{:}] = model.AutoDiffBackend.initVariablesAD(p, pm, wellVars{:});
    else
        wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
        % Dual porosity: adding matrix
        [p0, pm0, wellVars0{:}] = model.AutoDiffBackend.initVariablesAD(p0, pm0, wellVars0{:});  %#ok
    end
end
primaryVars = {'pressure', 'pressure_matrix', wellVarNames{:}};
gdz   = s.Grad(G.cells.centroids) * model.getGravityVector()';
%--------------------
%check for p-dependent tran mult:
trMult = 1;
if isfield(f, 'tranMultR'), trMult = f.tranMultR(p); end

%check for p-dependent porv mult:
pvMult = 1; pvMult0 = 1;
pvMultm = 1; pvMultm0 = 1;
if isfield(f, 'pvMultR')
    pvMult =  f.pvMultR(p);
    pvMult0 = f.pvMultR(p0);
end

if isfield(fm, 'pvMultR')
    pvMultm =  fm.pvMultR(pm);
    pvMultm0 = fm.pvMultR(pm0);
end
transMult=1;
if isfield(f, 'transMult')
   transMult=f.transMult(p); 
end

trans=s.T.*transMult;
% -------------------------------------------------------------------------
% water props (calculated at oil pressure OK?)
bW     = f.bW(p);
rhoW   = bW.*f.rhoWS;
% rhoW on face, avarge of neighboring cells (E100, not E300)
rhoWf  = s.faceAvg(rhoW);
mobW   = trMult./f.muW(p);
dpW     = s.Grad(p) - rhoWf.*gdz;
% water upstream-index
upcw = (value(dpW)<=0);
vW = - s.faceUpstr(upcw, mobW).*trans.*dpW;
bWvW = s.faceUpstr(upcw, bW).*vW;

% Matrix props
bWm = fm.bW(pm);

if model.outputFluxes
    state = model.storeFluxes(state, vW, [], []);
end

if model.extraStateOutput
    state = model.storebfactors(state, bW, [], []);
    state = model.storeMobilities(state, mobW, [], []);
    state = model.storeUpstreamIndices(state, upcw, [], []);
end
%% Transfer
vb = model.G.cells.volumes;

matrix_fields.pm = pm;
fracture_fields.pf = p;

transfer_model = model.transfer_model_object;

[Talpha] = transfer_model.calculate_transfer(model,fracture_fields,matrix_fields);

Tw = vb.*Talpha{1};

% EQUATIONS ---------------------------------------------------------------
names = {'water','water_matrix'};
types = {'cell','cell'};

% water:
water_fracture = (s.pv/dt).*( pvMult.*bW - pvMult0.*f.bW(p0) ) + s.Div(bWvW) + Tw;
water_matrix = (s.pv_matrix/dt).*( pvMultm.*bWm - pvMultm0.*fm.bW(pm0) ) - Tw;

eqs{1} = water_fracture;
eqs{2} = water_matrix;

% Dummy saturation
sW = ones(model.G.cells.num, 1);
[eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                                 {p}, {sW}, {mobW}, {rhoW}, ...
                                                                 {}, {}, ...
                                                                 drivingForces);

% well equations
[eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, types, wellSol0, wellSol, wellVars, wellMap,...
        p, {mobW}, {rhoW}, {}, {}, dt, opt);
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end
%--------------------------------------------------------------------------

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







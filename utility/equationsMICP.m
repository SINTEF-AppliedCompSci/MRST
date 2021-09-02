function [problem, state] = equationsMICP(state0, state, model, dt, ...
                                                   drivingForces, varargin)
% Assemble the linearized equations for the MICP system.
% 
% This function is modified from a file in The MATLAB Reservoir Simulation
% Toolbox (MRST), see
%   mrst/modules/ad-eor/utils/equationsOilWaterPolymer.m
%
% We refer to that function for a complete commented version of the file. 
% In this file we comment on the new added lines. 

%{
Partial copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
Partial copyright 2021 NORCE Norwegian Research Centre AS, Computational 
Geosciences and Modeling.

This file is part of the ad-micp module.

ad-micp is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ad-micp is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%}

opt = struct('Verbose', mrstVerbose, 'reverseMode', false, 'resOnly', ...
                                                   false, 'iteration', -1);        
opt = merge_options(opt, varargin{:});

% Model parameters
rhoW = model.fluid.rhoWS * model.fluid.cells;
mobW = model.fluid.cells / model.fluid.muw;
rho_b = model.fluid.rho_b;   
rho_c = model.fluid.rho_c; 
k_str = model.fluid.k_str;
mu_u = model.fluid.mu_u;                                                                                                               
k_a = model.fluid.k_a;    
k_d = model.fluid.k_d;
k_o = model.fluid.k_o;
k_u = model.fluid.k_u;
Yuc = model.fluid.Yuc;
mu = model.fluid.mu;  
Y = model.fluid.Y;
F = model.fluid.F;

% Properties at current timestep
[p, sW, m, o, u, b, c, wellSol] = model.getProps(state, 'pressure', ...
      'water', 'microorganism', 'oxygen', 'urea', 'biofilm', 'calcite', ...
                                                                'wellSol');
% Properties at previous timestep
[sW0, m0, o0, u0, b0, c0, wellSol0] = model.getProps(state0, 'water', ...
       'microorganism', 'oxygen', 'urea', 'biofilm', 'calcite', 'wellSol');

[wellVars, wellVarNames, wellMap] =...
                       model.FacilityModel.getAllPrimaryVariables(wellSol);
% Initialize independent variables.
if ~opt.resOnly
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode
        [p, sW, m, o, u, b, c, wellVars{:}] = ...
            initVariablesADI(p, sW, m, o, u, b, c, wellVars{:});
        primaryVars = [{'pressure'}, {'sW'}, {'microorganism'}, ...
         {'oxygen'}, {'urea'}, {'biofilm'}, {'calcite'}, wellVarNames(:)'];
    else
        [sW0, m0, o0, u0, b0, c0] = ...
            initVariablesADI(sW0, m0, o0, u0, b0, c0); 
        primaryVars = {'pressure', 'sW', 'microorganism', 'oxygen', ...
                                             'urea', 'biofilm', 'calcite'};
    end
else
    primaryVars = {'pressure', 'sW', 'microorganism', 'oxygen', ...
                                             'urea', 'biofilm', 'calcite'};
end
% Single-phase flow (sW = 1, sO = 0) (possible to extend to two-phase flow)
sO  = 1 - sW; 
gdz = model.getGravityGradient();

poro0 = model.rock.poro - c0 - b0;
poro = model.rock.poro - c - b;
if isobject(poro) == 1
        poro = poro.val;
end

[vW, vM, vO, vU, model] = getFluxAndPropsMICP(model, p, ...
                                                       m, o, u, gdz, poro);
state = model.storeFluxes(state, vW, 0 * vW, []);
[d_m, d_o, d_u, dpW] = getDispersionAnddpWMICP(model, state, poro);

op = model.operators;
pv = model.G.cells.volumes;
% Conservation of mass for water
water = (pv ./ dt) .* (poro .* sW - poro0 .* sW0) + op.Div(vW);
% Conservation of mass for CO2 (possible to extend to two-phase flow)
oil = sO;
% Conservation of mass for microorganisms
microorganism = (pv ./ dt) .* (poro .* m - poro0 .* m0) + ...
  op.Div(vM - d_m .* op.Grad(m)) - pv .* (m .* poro .* (Y * mu .* (o ./ ...
   (k_o + o)) - k_d - k_a) + b .* rho_b .* k_str .* (poro .* dpW) .^ 0.58);
% Conservation of mass for oxygen
oxygen = (pv ./ dt) .* (poro .* o - poro0 .* o0) + ...
op.Div(vO - d_o .* op.Grad(o)) + pv .* (m .* poro + rho_b .* b) .* F .* ...
                                                    mu .* (o ./ (k_o + o));
% Conservation of mass for urea
urea = (pv ./ dt) .* (poro .* u - poro0 .* u0) + ...
      op.Div(vU - d_u .* op.Grad(u)) + pv .* rho_b .* b * mu_u .* (u ./ ...
                                                                (k_u + u));
% Conservation of mass for biofilm
biofilm = (1 / dt) .* (b - b0) - ...
b .* (Y * mu .* o ./ (k_o + o) - k_d - rho_b .* b .* Yuc .* mu_u .* ...
                   (u ./ (k_u + u)) ./ (rho_c .* (poro + b)) - k_str .* ...
                         (poro .* dpW) .^ 0.58) - m .* poro .* k_a / rho_b;
% Conservation of mass for calcite
calcite = (1 / dt) .* (c - c0) - ...
                    rho_b .* b .* Yuc .* mu_u .* (u ./ (k_u + u)) ./ rho_c;                                                    

eqs   = {water, oil, microorganism, oxygen, urea, biofilm, calcite};
names = {'water', 'oil', 'microorganism', 'oxygen', 'urea', 'biofilm', ...
                                                                'calcite'};
types = {'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell'};

rho = {rhoW, rhoW};
mob = {mobW, 0 * mobW};
sat = {sW, sO};

[eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, ...
  types, state, {p, p}, sat, mob, rho, {}, {m, o, u, b, c}, drivingForces);

[eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, ...
   names, types, wellSol0, wellSol, wellVars, wellMap, p, mob, rho, {}, ...
                                                 {m, o, u, b, c}, dt, opt);                                           
                                            
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end
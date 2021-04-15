function [problem, state] = pressureEquationOilWaterPolymer(state0, state, model, dt, drivingForces, varargin)
%Undocumented Utility Function

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

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'propsPressure', [], ...
             'staticWells',  false, ...
             'iteration', -1);

opt = merge_options(opt, varargin{:});
assert(isempty(drivingForces.bc) && isempty(drivingForces.src))

W     = drivingForces.W;
op    = model.operators;
fluid = model.fluid;    

% Polymer shear thinning/thickening
usingShear = isfield(fluid, 'plyshearMult');

% Properties at current timestep
[p, sW, c, cmax, wellSol] = model.getProps(state, 'pressure', 'water', ...
    'polymer', 'polymermax', 'wellsol');

% Properties at previous timestep
[p0, sW0, c0, cmax0, wellSol0] = model.getProps(state0, 'pressure', 'water', ...
   'polymer', 'polymermax', 'wellsol');

%Initialization of independent variables ----------------------------------

[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);
if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [p, wellVars{:}] = ...
            model.AutoDiffBackend.initVariablesAD(p, wellVars{:});
    else
        assert(0, 'Backwards solver not supported for splitting');
    end
end
primaryVars = {'pressure', wellVarNames{:}};

p_prop = opt.propsPressure;
otherPropPressure = ~isempty(p_prop);
if ~otherPropPressure
    p_prop = p;
end

% -------------------------------------------------------------------------
sO  = 1 - sW;
sO0 = 1 - sW0;

[krW, krO] = model.evaluateRelPerm({sW, sO});

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, ...
    p_prop, p0);

% Modifiy relperm by mobility multiplier (if any)
krW = mobMult.*krW; krO = mobMult.*krO;

% Compute transmissibility
T = op.T.*transMult;

% Gravity contribution
gdz = model.getGravityGradient();

% Evaluate water and polymer properties
ads  = effads(c, cmax, model);
[vW, vP, bW, muWMult, mobW, mobP, rhoW, pW, upcw] = ...
    getFluxAndPropsWaterPolymer_BO(model, p_prop, sW, c, ads, ...
    krW, T, gdz); %, 'shear', ~otherPropPressure);
bW0 = fluid.bW(p0);
ads0  = effads(c0, cmax0, model);

% Evaluate oil properties
[vO, bO, mobO, rhoO, pO, upco, dpO] = getFluxAndPropsOil_BO(model, ...
    p_prop, sO, krO, T, gdz);
bO0 = getbO_BO(model, p0);

if otherPropPressure
    % We have used a different pressure for property evaluation, undo the
    % effects of this on the fluxes.
    dp_diff = op.Grad(p) - op.Grad(p_prop);
    
    vW = -op.faceUpstr(upcw, mobW).*op.T.*(dpW + dp_diff);
    vO = -op.faceUpstr(upco, mobO).*op.T.*(dpO + dp_diff);
    vP = -op.faceUpstr(upcw, mobP).*op.T.*(dpW + dp_diff);
    
    % If other property pressure is used, then the shear thinning is
    % computed here, and not above
    if usingShear
        poroFace  = op.faceAvg(model.rock.poro);
        faceArea  = model.G.faces.areas(op.internalConn);
        Vw        = vW./(poroFace .* faceArea); % water velocity
        muWMultf  = op.faceUpstr(upcw, muWMult);
        [shearMult, ~, shearReport] = ...
            getPolymerShearMultiplier(model, Vw, muWMultf);
        vW        = vW .* shearMult;
        vP        = vP .* shearMult;
        extraOutput.shearMult   = shearMult;
        extraOutput.shearReport = shearReport;
    end
end

% These are needed in transport solver, so we output them regardless of
% any flags set in the model.
state = model.storeFluxes(state, vW, vO, vP);
state = model.storeUpstreamIndices(state, upcw, upco, []);
if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, []);
    state = model.storeMobilities(state, mobW, mobO, mobP);
end

if usingShear
    % The shear multipliers are stored for use in the transport solver
    state = model.storeShearMultiplier(state, extraOutput.shearMult);
end


% EQUATIONS ---------------------------------------------------------------
% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.
bOvO = op.faceUpstr(upco, bO).*vO;
bWvW = op.faceUpstr(upcw, bW).*vW;
bWvP = op.faceUpstr(upcw, bW).*vP;

% Conservation of mass for water
water = (op.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + op.Div(bWvW);

% Conservation of mass for oil
oil = (op.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) + op.Div(bOvO);

% Conservation of polymer in water:
poro = model.rock.poro;
polymer = (op.pv.*(1-fluid.dps)/dt).*(pvMult.*bW.*sW.*c - ...
   pvMult0.*bW0.*sW0.*c0) + (op.pv/dt).* ...
   (fluid.rhoR.*((1-poro)./poro).*(ads-ads0) ) + op.Div(bWvP);

eqs   = {water, oil, polymer};
names = {'water', 'oil', 'polymer'};
types = {'cell', 'cell', 'cell'};

rho = {rhoW, rhoO};
mob = {mobW, mobO};
sat = {sW, sO};
[eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                                 {pW, p}, sat, mob, rho, ...
                                                                 {}, {c}, ...
                                                                 drivingForces);
                                                             
% Finally, add in and setup well equations
dissolved = {};
[eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, ...
    types, wellSol0, wellSol, wellVars, wellMap, p, mob, rho, dissolved, {c}, dt, opt);

eqs{1}   = (dt./op.pv).*(eqs{1}./bW + eqs{2}./bO);
names{1} = 'pressure';
types{1} = 'cell';

eqs   = eqs([1, 4:end]);
names = names([1, 4:end]);
types = types([1, 4:end]);

state.timestep = dt;
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

end

%--------------------------------------------------------------------------
% Helper functions
%--------------------------------------------------------------------------

% Effective adsorption, depending of desorption or not
function y = effads(c, cmax, model)
   if model.fluid.adsInx == 2
      y = model.fluid.ads(max(c, cmax));
   else
      y = model.fluid.ads(c);
   end
end

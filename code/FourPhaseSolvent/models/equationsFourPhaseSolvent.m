function [problem, state] = equationsFourPhaseSolvent(state0, state, model, dt, drivingForces, varargin)
% Equations for the four-phase solvent model

%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

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
             'iteration', -1);

opt = merge_options(opt, varargin{:});
opt.resOnly = false;

W = drivingForces.W;
op = model.operators;

% Properties at current timestep
[p, sW, sO, sG, wellSol] = model.getProps(state, 'pressure', 'water', ...
    'oil', 'gas', 'wellSol');

% Properties at previous timestep
[p0, sW0, sO0, sG0, wellSol0] = model.getProps(state0, 'pressure', 'water', ...
   'oil', 'gas', 'wellSol');

% Well variables
[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

if ~opt.resOnly
    if ~opt.reverseMode
        % define primary varible x and initialize
        [p, sW, sO, sG, wellVars{:}] = ...
            initVariablesADI(p, sW, sO, sG, wellVars{:});
    else
        % Set initial gradient to zero
        zw = zeros(size(bhp));
        [p0, sW0, sO0, sG0, zw, zw, zw, zw] = ...
            initVariablesADI(p0, sW0, sO0, sG0, zw, zw, zw, zw); %#ok
        clear zw;
    end
end

% We will solve for pressure, and water/oil/gas saturations (solvent
% saturation follows via the definition of saturations), and well rates +
% bhp.
primaryVars = {'pressure', 'sW', 'sO', 'sG', wellVarNames{:}};

sS  = 1 - sW  - sO  - sG ;
sS0 = 1 - sW0 - sO0 - sG0;

fluid = model.fluid;

% Get multipliers
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);
T = op.T.*transMult;

% Calculate residual saturations
[sWres, sOres , sSGres ] = computeResidualSaturations(fluid, p , sG , sS );
[~    , sO0res, sSG0res] = computeResidualSaturations(fluid, p0, sG0, sS0);

 % Calculate effective relperms
[krW, krO, krG, krS] = computeRelPermSolvent(fluid, p, sW, sO, sG, sS, sWres, sOres, sSGres, mobMult);

 % Calulate effective viscosities and densities
[muW, muO, muG, muS, rhoW, rhoO , rhoG , rhoS ] ...
    = computeViscositiesAndDensities(fluid, p , sO , sG , sS , sOres , sSGres );
[~  , ~  , ~  , ~  , ~   , rhoO0, rhoG0, rhoS0] ...
    = computeViscositiesAndDensities(fluid, p0, sO0, sG0, sS0, sO0res, sSG0res);

% Calulcate effective formation volume factors
[bW , bO , bG , bS ] = computeFormationVolumeFactors(fluid, p , rhoO , rhoG , rhoS );
[bW0, bO0, bG0, bS0] = computeFormationVolumeFactors(fluid, p0, rhoO0, rhoG0, rhoS0);

gdz = model.getGravityGradient();

[vW, vO, vG, vS, mobW, mobO, mobG, mobS, upcW, upcO, upcG, upcS] ...
    = getFluxAndPropsSolvent(fluid, p, krW, krO, krG, krS, muW, muO, muG, muS, rhoW, rhoO, rhoG, rhoS, T, gdz, op);

if model.outputFluxes
    state = model.storeFluxes(state, vW, vO, vG, vS);
end
if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, bG, bS);
    state = model.storeDensity(state, rhoW, rhoO, rhoG, rhoS);
    state = model.storeMobilities(state, mobW, mobO, mobG, mobS);
    state = model.storeUpstreamIndices(state, upcW, upcO, upcG, upcS);
end

% EQUATIONS ---------------------------------------------------------------
% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.
bWvW = op.faceUpstr(upcW, bW).*vW;
bOvO = op.faceUpstr(upcO, bO).*vO;
bGvG = op.faceUpstr(upcG, bG).*vG;
bSvS = op.faceUpstr(upcS, bS).*vS;

% Conservation of mass for water
water = (op.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + op.Div(bWvW);

% Conservation of mass for oil
oil = (op.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) + op.Div(bOvO);

% Conservation of mass for gas
gas = (op.pv/dt).*( pvMult.*bG.*sG - pvMult0.*bG0.*sG0 ) + op.Div(bGvG);

% Conservation of mass for solvent
solvent = (op.pv/dt).*( pvMult.*bS.*sS - pvMult0.*bS0.*sS0 ) + op.Div(bSvS);

if 0
    acc = zeros(model.G.cells.num,4);
    acc(:,1) = (op.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 );
    acc(:,2) = (op.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 );
    acc(:,3) = (op.pv/dt).*( pvMult.*bG.*sG - pvMult0.*bG0.*sG0 );
    acc(:,4) = (op.pv/dt).*( pvMult.*bS.*sS - pvMult0.*bS0.*sS0 );
    state.acc = acc;
    if isa(sOres, 'ADI')
        state.sOres = sOres.val;
    else
        state.sOres = sOres;
    end
end


eqs   = {water, oil, gas, solvent};
names = {'water', 'oil', 'gas', 'solvent'};
types = {'cell', 'cell', 'cell', 'cell'};

rho = {rhoW, rhoO, rhoG, rhoS};
mob = {mobW, mobO, mobG, mobS};
sat = {sW, sO, sG, sS};


pW = p;
[eqs, ~, qRes] = addFluxesFromSourcesAndBC(model, eqs, ...
                                       {pW, p},...
                                       rho, ...
                                       mob, ...
                                       sat, ...
                                       drivingForces);
                                   
if model.outputFluxes
    state = model.storeBoundaryFluxes(state, qRes{1}, qRes{2}, [], drivingForces);
end

% Finally, add in and setup well equations
if ~isempty(W)
    
    % Density of injected fluids are calculated using the 'saturations'
    % given in compi
    nc = arrayfun(@(w) numel(w.cells), W);
    sign = rldecode(vertcat(wellSol.sign), nc, 1);

    wc = vertcat(W.cells);
    wc = wc(sign>0);

    compi = rldecode(vertcat(W.compi), nc, 1);
    compi = compi(sign>0,:);

    muWell = cell(4,1);
    rhoWell = cell(4,1);
    [krWell{1}, krWell{2}, krWell{3}, krWell{4}] ...
        = computeRelPermSolvent(fluid, p(wc), compi(:,1), compi(:,2), compi(:,3), compi(:,4), 0,0,0, mobMult);
    [muWell{1}, muWell{2}, muWell{3}, muWell{4}, rhoWell{1}, rhoWell{2}, rhoWell{3}, rhoWell{4}] ...
        = computeViscositiesAndDensities(fluid, p(wc), compi(:,2), compi(:,3), compi(:,4), 0, 0);

    for i = 1:4
        rho{i}(wc) = rhoWell{i};
        mob{i}(wc) = mobMult*krWell{i}./muWell{i};
    end 
    
%     wm = model.FacilityModel;
    if ~opt.reverseMode
        [eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, types, wellSol0, wellSol, wellVars, wellMap, p, mob, rho, {}, {}, dt, opt);
    else
%         [eqs(4:7), names(4:7), types(4:7)] = wm.createReverseModeWellEquations(model, state0.wellSol, p0);
    end
end

% Scale solvent equations (this should be changed to CNV/MB convergence...)
eqs{4} = eqs{4}.*(dt./op.pv);

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

end
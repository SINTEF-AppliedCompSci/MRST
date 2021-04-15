function [problem, state] = equationsMiscibleOilWaterSolvent(state0, state, model, dt, ...
                                                     drivingForces, varargin)
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
             'iteration', -1);

opt = merge_options(opt, varargin{:});
opt.resOnly = false;

W = drivingForces.W;
s = model.operators;

% Properties at current timestep
[p, sW, sO, wellSol] = model.getProps(state, 'pressure', 'water', ...
                                                         'oil', 'wellSol');

% Properties at previous timestep
[p0, sW0, sO0, wellSol0] = model.getProps(state0, 'pressure', 'water', ...
                                                         'oil', 'wellSol');

[qWell, bhp, wellVars, wellVarNames, wellMap] = ...
                       model.FacilityModel.getAllPrimaryVariables(wellSol);

if ~opt.resOnly
    if ~opt.reverseMode
        % define primary varible x and initialize
        [p, sW, sO, qWell{:}, bhp] = ...
            initVariablesADI(p, sW, sO, qWell{:}, bhp);
    else
        % Set initial gradient to zero
        zw = zeros(size(bhp));
        [p0, sW0, sO0, zw, zw, zw, zw] = ...
            initVariablesADI(p0, sW0, sO0, zw, zw, zw, zw); %#ok
        clear zw;
    end
end

% We will solve for pressure, water saturation and oil saturation (solvent
% saturation follows via the definition of saturations),  and well rates +
% bhp.
primaryVars = {'pressure', 'sW', 'sO', wellVarNames{:}};

sG  = 1 - sW  - sO ;
sG0 = 1 - sW0 - sO0;

% Get dynamic quantities
[kr, mu, rho, b, b0, pvMult, pvMult0, T] ...
    = getDynamicQuantitiesMiscibleOilWaterSolvent(model, p0, p, sW, sO, sG, sO0, sG0);

krW  = kr{1} ; krO  = kr{2} ; krG  = kr{3} ;
rhoW = rho{1}; rhoO = rho{2}; rhoG = rho{3};
muW  = mu{1} ; muO  = mu{2} ; muG  = mu{3} ;
bW   = b{1}  ; bO   = b{2}  ; bG   = b{3}  ;
bW0  = b0{1} ; bO0  = b0{2} ; bG0  = b0{3} ;

gdz = model.getGravityGradient();
op = model.operators;

[vW, mobW, upcW] = getFlux(p, rhoW, krW, muW, T, gdz, op);
[vO, mobO, upcO] = getFlux(p, rhoO, krO, muO, T, gdz, op);
[vG, mobG, upcG] = getFlux(p, rhoG, krG, muG, T, gdz, op);


if model.outputFluxes
    state = model.storeFluxes(state, vW, vO, vG);
end
if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, bG, bS);
    state = model.storeMobilities(state, mobW, mobO, mobG, mobS);
    state = model.storeUpstreamIndices(state, upcW, upcO, upcG, upcS);
end

% EQUATIONS ---------------------------------------------------------------
% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.
bWvW = s.faceUpstr(upcW, bW).*vW;
bOvO = s.faceUpstr(upcO, bO).*vO;
bGvG = s.faceUpstr(upcG, bG).*vG;

% Conservation of mass for water
water = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + s.Div(bWvW);

% Conservation of mass for oil
oil = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) + s.Div(bOvO);

% Conservation of mass for gas
gas = (s.pv/dt).*( pvMult.*bG.*sG - pvMult0.*bG0.*sG0 ) + s.Div(bGvG);

eqs   = {water, oil, gas};
names = {'water', 'oil', 'gas'};
types = {'cell', 'cell', 'cell'};

% Add in any fluxes / source terms prescribed as boundary conditions.

rho = {rhoW, rhoO, rhoG};
mob = {mobW, mobO, mobG};
sat = {sW, sO, sG};

wc = W.cells;
wcInj = wc([wellSol.sign] == 1);
compi = reshape([W.compi],3, [])';
compi = compi([wellSol.sign] == 1,:);

% Get dynamic quantities
[~, ~, rhoWell, ~, ~, ~, ~, ~] ...
    = getDynamicQuantitiesMiscibleOilWaterSolvent(model, 0, p(wcInj), compi(:,1), compi(:,2), compi(:,3), 0,0);

for i = 1:3
    rho{i}.val(wcInj) = rhoWell{i}.val;
end

% [eqs, ~, qRes] = addFluxesFromSourcesAndBC(model, eqs, ...
%                                        {pW, p},...
%                                        rho, ...
%                                        mob, ...
%                                        sat, ...
%                                        drivingForces);
% if model.outputFluxes
%     state = model.storeBoundaryFluxes(state, qRes{1}, qRes{2}, [], drivingForces);
% end
% Finally, add in and setup well equations

if ~isempty(W)
    wm = model.FacilityModel;
    if ~opt.reverseMode
        [eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, ...
                                                          names, types, wellSol0, ...
                                                          wellSol, qWell, bhp, ...
                                                          wellVars, wellMap, ...
                                                          p, mob, rho, {}, {}, ...
                                                          dt, opt);
    else
        [eqs(4:7), names(4:7), types(4:7)] = wm.createReverseModeWellEquations(model, ...
                                                          state0.wellSol, p0);
    end
end
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end


%--------------------------------------------------------------------------

function [v, mob, upc] = getFlux(p, rho, kr, mu, T, gdz, op)

    rhof  = op.faceAvg(rho);
    mob   = kr./mu;
    dp    = op.Grad(p) - rhof.*gdz;
    
    upc  = (double(dp)<=0);
    v   = - op.faceUpstr(upc, mob).*T.*dp;
    
end

function [problem, state] = equationsOilWaterSolvent(state0, state, model, dt, ...
                                                     drivingForces, varargin)
                                                 
                                                 
opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'iteration', -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.W;
s = model.operators;

% Properties at current timestep
[p, sW, sO, wellSol] = model.getProps(state, 'pressure', 'water', ...
    'oil', 'wellSol');

% Properties at previous timestep
[p0, sW0, sO0, wellSol0] = model.getProps(state0, 'pressure', 'water', ...
   'oil', 'wellSol');

[qWell, pBH, wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);
% Initialize independent variables.
if ~opt.resOnly
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode
        [p, sW, sO, qWell{:}, pBH, wellVars{:}] = ...
            initVariablesADI(p, sW, sO, qWell{:}, pBH, wellVars{:});
    else
        zw = zeros(size(pBH));
        [p0, sW0, sO0, zw, zw, zw, zw] = ...
            initVariablesADI(p0, sW0, sO0, zw, zw, zw, zw); %#ok
        clear zw
    end
end

% We will solve for pressure, water saturation and oil saturation (solvent
% saturation follows via the definition of saturations),  and well rates +
% bhp.
primaryVars = {'pressure', 'sW', 'sO', wellVarNames{:}};

sS  = 1 - sW - sO;
sS0 = 1 - sW0 - sO0;

% Get dynamic quantities
[krW , krO , krG , krS , ...
 muW , muO , muG , muS , ...
 rhoW, rhoO, rhoG, rhoS, ...
 bW  , bO  , bG  , bS  , ...
 bW0 , bO0 , bG0 , bS0 , ...
 pvMult, transMult, mobMult, pvMult0, T] ...
               = getDynamicQuantitiesSolvent(model, p0, p, sW, sO, [], sS);

gdz = model.getGravityGradient();
op = model.operators;
[vW, mobW, upcW] = getFlux(p, rhoW, krW, muW, T, gdz, op);
[vO, mobO, upcO] = getFlux(p, rhoO, krO, muO, T, gdz, op);
[vS, mobS, upcS] = getFlux(p, rhoS, krS, muS, T, gdz, op);


if model.outputFluxes
    state = model.storeFluxes(state, vW, vO, vS);
end
if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, bS);
    state = model.storeMobilities(state, mobW, mobO, mobS);
    state = model.storeUpstreamIndices(state, upcW, upcO, upcS);
end

% EQUATIONS ---------------------------------------------------------------
% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.
bWvW = s.faceUpstr(upcW, bW).*vW;
bOvO = s.faceUpstr(upcO, bO).*vO;
bSvS = s.faceUpstr(upcS, bS).*vS;

% Conservation of mass for water
water = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + s.Div(bWvW);

% Conservation of mass for oil
oil = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) + s.Div(bOvO);

% Conservation of mass for solvent
solvent = (s.pv/dt).*( pvMult.*bS.*sS - pvMult0.*bS0.*sS0 ) + s.Div(bSvS);

eqs   = {water, oil, solvent};
names = {'water', 'oil', 'solvent'};
types = {'cell', 'cell', 'cell'};

% % Add in any fluxes / source terms prescribed as boundary conditions.
rho = {rhoW, rhoO, rhoS};
mob = {mobW, mobO, mobS};
sat = {sW, sO, sS};
% [eqs, qBC, qRes, BCTocellMap, qSRC, srcCells] = addFluxesFromSourcesAndBC(...
%    model, eqs, {pW, p}, rho, mob,  ...
%    sat, drivingForces);

% if model.outputFluxes
%     state = model.storeBoundaryFluxes(state, qRes{1}, qRes{2}, [], drivingForces);
% end
% Add polymer boundary conditions

% Finally, add in and setup well equations
if ~isempty(W) 
    wm = model.FacilityModel;
    if ~opt.reverseMode
        [eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, types, wellSol0, wellSol, qWell, pBH, wellVars, wellMap, p, mob, rho, {}, {}, dt, opt);
    else
        [eqs(3:5), names(3:5), types(3:5)] = wm.createReverseModeWellEquations(model, state0.wellSol, p0);
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
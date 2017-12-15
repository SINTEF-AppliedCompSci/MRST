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
                                                 
opt = struct('Verbose'    , mrstVerbose, ...
             'reverseMode', false      ,...
             'resOnly'    , false      ,...
             'iteration'  , -1         );

opt = merge_options(opt, varargin{:});
opt.resOnly = false;

W     = drivingForces.W;
op    = model.operators;
fluid = model.fluid;

% Properties at current timestep
[p , sW , sG , sS , rs , rv , wellSol ] = model.getProps(state, ...
             'pressure', 'water', 'gas', 'solvent', 'rs', 'rv', 'wellSol');
% Properties at previous timestep
[p0, sW0, sG0, sS0, rs0, rv0, wellSol0] = model.getProps(state0, ...
             'pressure', 'water', 'gas', 'solvent', 'rs', 'rv', 'wellSol');

% Well variables
[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

solventInWater = true;

% Initialize primary variables (gas saturation taken to be sG + sS)
if solventInWater
    st  = model.getCellStatusVO(state,  1-(sW + sG  + sS) , sW  + sS , sG );
    st0 = model.getCellStatusVO(state0, 1-(sW0+ sG0 + sS0), sW0 + sS0, sG0);
else
    st  = model.getCellStatusVO(state,  1-(sW + sG  + sS) , sW , sG  + sS );
    st0 = model.getCellStatusVO(state0, 1-(sW0+ sG0 + sS0), sW0, sG0 + sS0);
end

if model.disgas || model.vapoil
    % X is either Rs, Rv or Sg, depending on each cell's saturation status
    x = st{1}.*rs + st{2}.*rv + st{3}.*sG;
    gvar = 'x';
else
    x = sG;
    gvar = 'sG';
end

if ~opt.resOnly,
    if ~opt.reverseMode,
        % define primary varible x and initialize
        [p, sW, x, sS, wellVars{:}] = initVariablesADI(p, sW, x, sS, wellVars{:});
    else
        wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
        x0 = st0{1}.*rs0 + st0{2}.*rv0 + st0{3}.*sG0;
        [p0, sW0, x0, sS0, wellVars0{:}] = initVariablesADI(p0, sW0, x0, sS0, wellVars0{:}); %#ok
        if solventInWater
            [sG0, rs0, rv0] = calculateHydrocarbonsFromStatusBO(model, st0, 1-(sW0 + sS0), x0, rs0, rv0, p0);
        else
            [sG0, rs0, rv0] = calculateHydrocarbonsFromStatusBO(model, st0, 1- sW0, x0 + sS0, rs0, rv0, p0);
        end
    end
end
    
if ~opt.reverseMode
    % Compute values from status flags. If we are in reverse mode, these
    % values have already converged in the forward simulation.
    if solventInWater
        [sG, rs, rv, rsSat, rvSat] = calculateHydrocarbonsFromStatusBO(model, st, 1 - (sW + sS), x, rs, rv, p);
    else
        [sG, rs, rv, rsSat, rvSat] = calculateHydrocarbonsFromStatusBO(model, st, 1 - sW , x + sS, rs, rv, p);
    end
end

% We will solve for pressure, and water/gas/solvent saturations (oil
% saturation follows via the definition of saturations), and well rates +
% bhp.
primaryVars = {'pressure', 'sW', gvar, 'sS', wellVarNames{:}};

sO  = 1 - (sW  + sG  + sS );
sO0 = 1 - (sW0 + sG0 + sS0);

% Calculate residual saturations
[sWr, sOr , sGc ] = computeResidualSaturations(model, p , sW , sG , sS , state0);
if isfield(state0, 'sOr') && isfield(state0, 'sGcr')
    sOr0 = state0.sOr;
    sGc0 = state0.sGc;
else
    [~    , sOr0, sGc0] = computeResidualSaturations(model, p0, sW0, sG0, sS0, state0);
end

% Get multipliers
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);
T = op.T.*transMult;

 % Calculate effective relperms
[krW, krO, krG, krS] = computeRelPermSolvent(model, p, sW, sO, sG, sS, sWr, sOr, sGc, mobMult);

 % Calulate effective viscosities and densities
[muW, muO, muG, muS, rhoW, rhoO , rhoG , rhoS , bW , bO , bG , bS , pW, pG] ...
    = computeViscositiesAndDensities(model, p , sW, sO , sG , sS , sOr , sGc , rs, rv, ~st{1}, ~st{2});
[~  , ~  , ~  , ~  , ~   , rhoO0, rhoG0, rhoS0, bW0, bO0, bG0, bS0, ~ , ~] ...
    = computeViscositiesAndDensities(model, p0, sW0, sO0, sG0, sS0, sOr0, sGc0, rs0, rv0, ~st0{1}, ~st0{2});

% % Calulcate effective formation volume factors
% [bW , bO , bG , bS ] = computeFormationVolumeFactors(fluid, p , rhoO , rhoG , rhoS );
% [bW0, bO0, bG0, bS0] = computeFormationVolumeFactors(fluid, p0, rhoO0, rhoG0, rhoS0);

gdz = model.getGravityGradient();

[vW, vO, vG, vS, mobW, mobO, mobG, mobS, upcW, upcO, upcG, upcS] ...
    = getFluxAndPropsSolvent(fluid, p, krW, krO, krG, krS, muW, muO, muG, muS, rhoW, rhoO, rhoG, rhoS, sW, sG, sS, T, gdz, op);

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
if model.vapoil
    % The model allows oil to vaporize into the gas phase. The conservation
    % equation for oil must then include the fraction present in the gas
    % phase.
    rvbGvG = op.faceUpstr(upcG, rv).*bGvG;
    % Final equation
    oil = (op.pv/dt).*( pvMult .*(bO.* sO  + rv .*bG .*sG )-    ...
                        pvMult0.*(bO0.*sO0 + rv0.*bG0.*sG0) ) + ...
           op.Div(bOvO + rvbGvG);
else
    oil = (op.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) + op.Div(bOvO);
end

% Conservation of mass for gas
if model.disgas
    % The gas transported in the oil phase.
    rsbOvO = op.faceUpstr(upcO, rs).*bOvO;
    
    gas = (op.pv/dt).*( pvMult.* (bG.* sG  + rs.* bO.* sO) - ...
        pvMult0.*(bG0.*sG0 + rs0.*bO0.*sO0 ) ) + ...
        op.Div(bGvG + rsbOvO);
else
    gas = (op.pv/dt).*( pvMult.*bG.*sG - pvMult0.*bG0.*sG0 ) + op.Div(bGvG);
end

% Conservation of mass for solvent
solvent = (op.pv/dt).*( pvMult.*bS.*sS - pvMult0.*bS0.*sS0 ) + op.Div(bSvS);

if 1
    
    acc = zeros(model.G.cells.num,4);
    acc(:,1)    = (op.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 );
    acc(:,2)    = (op.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 );
    acc(:,3)    = (op.pv/dt).*( pvMult.*bG.*sG - pvMult0.*bG0.*sG0 );
    acc(:,4)    = (op.pv/dt).*( pvMult.*bS.*sS - pvMult0.*bS0.*sS0 );
    state.acc   = acc;
    
    flux = zeros(model.G.cells.num,4);
    flux(:,1)    = op.Div(bWvW);
    flux(:,2)    = op.Div(bOvO);
    flux(:,3)    = op.Div(bGvG);
    flux(:,4)    = op.Div(bSvS);
    state.flux   = flux;
    
end

state.sr = [double(sOr), double(sGc)];
state.sOr = double(sOr);
state.sGc = double(sGc);


eqs   = {water, oil, gas, solvent};
names = {'water', 'oil', 'gas', 'solvent'};
types = {'cell', 'cell', 'cell', 'cell'};

rho = {rhoW, rhoO, rhoG, rhoS};
mob = {mobW, mobO, mobG, mobS};
sat = {sW, sO, sG, sS};
dissolved = model.getDissolutionMatrix(rs, rv);
% 
% [eqs, ~, qRes] = addFluxesFromSourcesAndBC(model, eqs, ...
%                                        {pW, p, pG, pG},...
%                                        rho, ...
%                                        mob, ...
%                                        sat, ...
%                                        drivingForces);
%                                    
% if model.outputFluxes
%     state = model.storeBoundaryFluxes(state, qRes{1}, qRes{2}, [], drivingForces);
% end

% Finally, add in and setup well equations
if ~isempty(W)
   
if 0
    % Density of injected fluids are calculated using the 'saturations'
    % given in compi
    nc = arrayfun(@(w) numel(w.cells), W);
    sign = rldecode(vertcat(wellSol.sign), nc, 1);

    wc = vertcat(W.cells);
    wc = wc(sign>0);

    compi = rldecode(vertcat(W.compi), nc, 1);
    compi = compi(sign>0,:);

    [mobMultw, rsw, rvw, isSatOw, isSatGw] = getWellValue(wc, mobMult, rs, rv, ~st{1}, ~st{2}); 

    [muWell, rhoWell] = deal(cell(4,1));
    [krWell{1}, krWell{2}, krWell{3}, krWell{4}] ...
        = computeRelPermSolvent(model, p(wc), compi(:,1), compi(:,2), compi(:,3), compi(:,4), sWr, sOr(wc), sGc(wc), mobMultw);
    [muWell{1}, muWell{2}, muWell{3}, muWell{4}, rhoWell{1}, rhoWell{2}, rhoWell{3}, rhoWell{4}] ...
        = computeViscositiesAndDensities(model, p(wc), compi(:,1), compi(:,2), compi(:,3), compi(:,4), sOr(wc), sGc(wc), rsw, rvw, isSatOw, isSatGw);
    for nPh = 1:4
        rho{nPh}(wc) = rhoWell{nPh};
        mob{nPh}(wc) = krWell{nPh}./muWell{nPh};
    end 
    
elseif 0
    
    nc = arrayfun(@(w) numel(w.cells), W);
    sign = rldecode(vertcat(wellSol.sign), nc, 1);

    wc = vertcat(W.cells);
    wc = wc(sign>0);

    compi = rldecode(vertcat(W.compi), nc, 1);
    compi = compi(sign>0,:);

    [mobMultw, rsw, rvw, isSatOw, isSatGw] = getWellValue(wc, mobMult, rs, rv, ~st{1}, ~st{2}); 
    
    FSGw = fluid.satFrac(compi(:,4), compi(:,3) + compi(:,4));
    krGTw = fluid.krG(compi(:,3) + compi(:,4));
    krGw  = fluid.krFG(1-FSGw).*krGTw;
    krSw  = fluid.krFS(FSGw).*krGTw;
    
    krWell   = {mobMultw.*fluid.krW(compi(:,1)), mobMultw.*fluid.krO(compi(:,2)), ...
               mobMult.*krGw, mobMult.*krSw};
    fld = fluid; fld.mixPar = 0; fld.mixParRho = 0;
    mdl = model; mdl.fluid = fld;

    [muWell, rhoWell] = deal(cell(4,1));
    [muWell{1}, muWell{2}, muWell{3}, muWell{4}, rhoWell{1}, rhoWell{2}, rhoWell{3}, rhoWell{4}] ...
        = computeViscositiesAndDensities(mdl, p(wc), compi(:,1), compi(:,2), compi(:,3), compi(:,4), sOr(wc), sGc(wc), rsw, rvw, isSatOw, isSatGw);
    for nPh = 1:4
        rho{nPh}(wc) = rhoWell{nPh};
        mob{nPh}(wc) = krWell{nPh}./muWell{nPh};
    end 
    
elseif 0
    
    nc = arrayfun(@(w) numel(w.cells), W);
    sign = rldecode(vertcat(wellSol.sign), nc, 1);
    wc = vertcat(W.cells);
    wc = wc(sign>0);
    
    Mpw = fluid.Mp(p(wc));
    rhoS = model.getSurfaceDensities();
    
    bO_i = fluid.bO(p, rs, ~st{1});
    bG_i = fluid.bG(p);
    bS_i = fluid.bS(p);
    b = {bW, bO_i, bG_i, bS_i};
        
    for phNo = 1:4
        rho{phNo}(wc) = rho{phNo}(wc).*Mpw + rhoS(phNo).*b{phNo}(wc).*(1-Mpw);
    end
    
end
    
    [eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, types, wellSol0, wellSol, wellVars, wellMap, p, mob, rho, dissolved, {}, dt, opt);

end

% Scale solvent equations (this should be changed to CNV/MB convergence...)
eqs{4} = eqs{4}.*(dt./op.pv);

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

end

function varargout = getWellValue(wellCells, varargin)
    
    for i = 1:numel(varargin)
        v = varargin{i};
        if numel(double(v)) == 1
            varargout{i} = v;
        else
            varargout{i} = v(wellCells);
        end 
    end

end

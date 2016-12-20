function [problem, state] = equationsOilWaterPolymer(state0, state, model, dt, ...
                                                     drivingForces, varargin)
%
% SYNOPSIS:
%   function [problem, state] = equationsOilWaterPolymer(state0, state, model, dt, drivingForces, varargin)
%
% DESCRIPTION: Assemble the linearized equations for an oil-water-polymer
% system, computing both the residuals and the Jacobians. Returns the result as
% an instance of the class LinearizedProblem which can be solved using instances
% of LinearSolverAD.
%
% A description of the modeling equations can be found in the directory
% ad-eor/docs.
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
% SEE ALSO: LinearizedProblem, LinearSolverAD, equationsOilWater, OilWaterPolymerModel
%

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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


% Get linearized problem for oil/water/polymer system with black oil-style
% properties
opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'iteration', -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.W;
s = model.operators;
f = model.fluid;

% Shear thinning/thickening
usingShear = isfield(f, 'plyshearMult');

% Properties at current timestep
[p, sW, c, cmax, wellSol] = model.getProps(state, 'pressure', 'water', ...
    'polymer', 'polymermax', 'wellsol');

% Properties at previous timestep
[p0, sW0, c0, cmax0] = model.getProps(state0, 'pressure', 'water', ...
   'polymer', 'polymermax');

pBH    = vertcat(wellSol.bhp);
qWs    = vertcat(wellSol.qWs);
qOs    = vertcat(wellSol.qOs);
qWPoly = vertcat(wellSol.qWPoly);

% Initialize independent variables.
if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [p, sW, c, qWs, qOs, qWPoly, pBH] = ...
            initVariablesADI(p, sW, c, qWs, qOs, qWPoly, pBH);
    else
        zw = zeros(size(pBH));
        [p0, sW0, c0, zw, zw, zw, zw] = ...
            initVariablesADI(p0, sW0, c0, zw, zw, zw, zw); %#ok
        clear zw
    end
end

% We will solve for pressure, water saturation (oil saturation follows via
% the definition of saturations), polymer concentration and well rates +
% bhp.
primaryVars = {'pressure', 'sW', 'polymer', 'qWs', 'qOs', 'qWPoly', 'bhp'};

% Evaluate relative permeability
sO  = 1 - sW;
sO0 = 1 - sW0;

[krW, krO] = model.evaluateRelPerm({sW, sO});

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

% Modifiy relperm by mobility multiplier (if any)
krW = mobMult.*krW; krO = mobMult.*krO;

% Compute transmissibility
T = s.T.*transMult;

% Gravity contribution
gdz = model.getGravityGradient();

% Evaluate water and polymer properties
ads  = effads(c, cmax, model);
ads0 = effads(c0, cmax0, model);
[vW, vP, bW, muWMult, mobW, mobP, rhoW, pW, upcw] = ...
    getFluxAndPropsWaterPolymer_BO(model, p, sW, c, ads, ...
    krW, T, gdz);
bW0 = model.fluid.bW(p0);

% Evaluate oil properties
[vO, bO, mobO, rhoO, p, upco] = getFluxAndPropsOil_BO(model, p, sO, krO, T, gdz);
bO0 = getbO_BO(model, p0);

if model.outputFluxes
    state = model.storeFluxes(state, vW, vO, vP);
end

if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, []);
    state = model.storeMobilities(state, mobW, mobO, mobP);
    state = model.storeUpstreamIndices(state, upcw, upco, []);
end

if model.extraPolymerOutput
    state = model.storeEffectiveWaterVisc(state, extraOutput.muWeff);
    state = model.storeEffectivePolymerVisc(state, extraOutput.muPeff);
    state = model.storePolymerAdsorption(state, ads);
    state = model.storeRelpermReductionFactor(state, extraOutput.Rk);
    if usingShear
        state = model.storeShearMultiplier(state, extraOutput.shearMult);
        state.ShearThinningReport = extraOutput.shearReport;
    end
end


% EQUATIONS ---------------------------------------------------------------
% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.
bOvO = s.faceUpstr(upco, bO).*vO;
bWvW = s.faceUpstr(upcw, bW).*vW;
bWvP = s.faceUpstr(upcw, bW).*vP;

% Conservation of mass for water
water = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + s.Div(bWvW);

% Conservation of mass for oil
oil = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) + s.Div(bOvO);

% Conservation of polymer in water:
poro = model.rock.poro;
f    = model.fluid;
polymer = (s.pv.*(1-f.dps)/dt).*(pvMult.*bW.*sW.*c - ...
   pvMult0.*bW0.*sW0.*c0) + (s.pv/dt).* ...
   ( f.rhoR.*((1-poro)./poro).*(ads-ads0) ) + s.Div(bWvP);

eqs   = {water, oil, polymer};
names = {'water', 'oil', 'polymer'};
types = {'cell', 'cell', 'cell'};

% TODO: % Fix for (almost) zero water in the well
% if isa(poly, 'ADI')
%    epsilon = 1.e-8;
%    epsilon = sqrt(epsilon)*mean(abs(diag(poly.jac{2})));
%    bad     = abs(diag(poly.jac{2})) < epsilon;
%    poly(bad) = c(bad);
% end


% Add in any fluxes / source terms prescribed as boundary conditions.
[eqs, qBC, ~, BCTocellMap, qSRC, srcCells] = addFluxesFromSourcesAndBC(...
   model, eqs, {pW, p}, {rhoW, rhoO}, {mobW, mobO}, {bW, bO},  ...
   {sW, sO}, drivingForces);

% Add polymer boundary conditions
if ~isempty(drivingForces.bc) && isfield(drivingForces.bc, 'poly')
   injInx  = qBC{1} > 0; % water inflow indecies
   cbc     = (BCTocellMap')*c; % BCTocellMap' = cellToBCMap
   cbc(injInx) = drivingForces.bc.poly(injInx);
   eqs{3}  = eqs{3} - BCTocellMap*(cbc.*qBC{1});
end

% Add polymer source
if ~isempty(drivingForces.src) && isfield(drivingForces.src, 'poly')
   injInx  = qSRC{1} > 0; % water inflow indecies
   csrc    = c(srcCells);
   csrc(injInx) = drivingForces.src.poly(injInx);
   eqs{3}(srcCells) = eqs{3}(srcCells) - csrc.*qSRC{1};
end

% Finally, add in and setup well equations
if ~isempty(W)
    wm = model.wellmodel;
    if ~opt.reverseMode
        wc   = vertcat(W.cells);
        pw   = p(wc);
        rhos = [f.rhoWS, f.rhoOS];
        bw   = {bW(wc), bO(wc)};
        mw   = {mobW(wc), mobO(wc)};
        s    = {sW(wc), sO(wc)};
        
        % Polymer well equations
        [~, wciPoly, iInxW] = getWellPolymer(W);        
        if usingShear
            % Compute shear rate multiplier for wells
            % The water velocity is computed using a the reprensentative 
            % radius rR.
            % rR = sqrt(rW * rA)
            % rW is the well bore radius.
            % rA is the equivalent radius of the grid block in which the 
            %    wellis completed.
            
            assert(isfield(W, 'rR'), ...
                'The representative radius needs to be suppplied.');
            
            muWMultW = muWMult(wc);
            % Maybe should also apply this for PRODUCTION wells.
            muWMultW((iInxW(wciPoly==0))) = 1;
            
            cqs = wm.computeWellFlux(model, W, wellSol, pBH, ...
                {qWs, qOs}, pw, rhos, bw, mw, s, {},...
                'nonlinearIteration', opt.iteration);
            
            % The following formulations assume that the wells are always
            % in the z direction 
            % IMPROVED HERE LATER
            [~, ~, dz] = cellDims(model.G, wc);
            
            % HACK FOR 2D MODELS
            if all(dz==0)
                dz(:) = 1;
            end
            
            if model.extraPolymerOutput
                cqsW0 = double(cqs{1});
                mobW0 = double(mw{1});
            end
            
            rR = vertcat(W.rR);
            A  = rR.*dz*2*pi; % representative area of each well cell
            VW0W = double(bW(wc)).*double(cqs{1})./(poro(wc).*A);
            [shearMultW, VW1W] = getPolymerShearMultiplier(model, ...
                VW0W, muWMultW);
            
            % Apply shear velocity multiplier
            mw{1} = mw{1}.*shearMultW;
        end
        
        
        [cqs, weqs, ctrleqs, wc, state.wellSol] = ...
            wm.computeWellFlux(model, W, wellSol, ...
            pBH, {qWs, qOs}, pw, rhos, bw, mw, s, {},...
            'nonlinearIteration', opt.iteration);

        % Store the well equations (relate well bottom hole pressures to
        % influx).
        eqs(4:5) = weqs;
        % Store the control equations (trivial equations ensuring that each
        % well will have values corresponding to the prescribed value)
        eqs{7} = ctrleqs;
        % Add source terms to the equations. Negative sign may be
        % surprising if one is used to source terms on the right hand side,
        % but this is the equations on residual form.
        eqs{1}(wc) = eqs{1}(wc) - cqs{1};
        eqs{2}(wc) = eqs{2}(wc) - cqs{2};

        % Polymer well equations
        [~, wciPoly, iInxW] = getWellPolymer(W);
        cw        = c(wc);
        cw(iInxW) = wciPoly;
        cbarw     = cw/f.cmax;

        % Divide away water mobility and add in polymer
        %bWqP = cw.*cqs{1}./(a + (1-a).*cbarw);
        bWqP = cw.*cqs{1};
        eqs{3}(wc) = eqs{3}(wc) - bWqP;

        % Well polymer rate for each well is water rate in each perforation
        % multiplied with polymer concentration in that perforated cell.
        perf2well = getPerforationToWellMapping(W);
        Rw = sparse(perf2well, (1:numel(perf2well))', 1, ...
           numel(W), numel(perf2well));
        cqsPoly = Rw*(cqs{1}.*cw);
        eqs{6}  = qWPoly - cqsPoly;
        
        % Save extra polymer welldata if requested
        if model.extraPolymerOutput
            cqsPoly    = double(cqsPoly);
            if usingShear
                shearMultW = double(shearMultW);
            end
            for wnr = 1:numel(state.wellSol)
                ix = perf2well == wnr;
                state.wellSol(wnr).cqsPoly = cqsPoly(wnr);
                if usingShear
                    state.wellSol(wnr).shearMult = shearMultW(ix);
                end
            end
        end
        
        names(4:7) = {'waterWells', 'oilWells', 'polymerWells', ...
            'closureWells'};
        types(4:7) = {'perf', 'perf', 'perf', 'well'};
    else
        [eq, n, typ] = ...
            wm.createReverseModeWellEquations(model, state0.wellSol, p0);
        % Add another equation for polymer well rates
        [eqs{4:7}] = deal(eq{1});
        [names{4:7}] = deal(n{1});
        [types{4:7}] = deal(typ{1});
    end
end
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end


%--------------------------------------------------------------------------

% Effective adsorption, depending of desorption or not
% adsInx=1: The polymer adsorption isotherm is retraced whenever the local 
%           polymer concentration in the solution decreases.
% adsInx=2: No polymer desorption may occur.
function y = effads(c, cmax, model)
   if model.fluid.adsInx == 2
      y = model.fluid.ads(max(c, cmax));
   else
      y = model.fluid.ads(c);
   end
end


function [dx, dy, dz] = cellDims(G, ix)
% cellDims -- Compute physical dimensions of all cells in single well
%
% SYNOPSIS:
%   [dx, dy, dz] = cellDims(G, ix)
%
% PARAMETERS:
%   G  - Grid data structure.
%   ix - Cells for which to compute the physical dimensions
%
% RETURNS:
%   dx, dy, dz -- [dx(k) dy(k)] is bounding box in xy-plane, while dz(k) =
%                 V(k)/dx(k)*dy(k)

    n = numel(ix);
    [dx, dy, dz] = deal(zeros([n, 1]));

    ixc = G.cells.facePos;
    ixf = G.faces.nodePos;

    for k = 1 : n,
       c = ix(k);                                     % Current cell
       f = G.cells.faces(ixc(c) : ixc(c + 1) - 1, 1); % Faces on cell
       e = mcolon(ixf(f), ixf(f + 1) - 1);            % Edges on cell

       nodes  = unique(G.faces.nodes(e, 1));          % Unique nodes...
       coords = G.nodes.coords(nodes,:);            % ... and coordinates

       % Compute bounding box
       m = min(coords);
       M = max(coords);

       % Size of bounding box
       dx(k) = M(1) - m(1);
       if size(G.nodes.coords, 2) > 1,
          dy(k) = M(2) - m(2);
       else
          dy(k) = 1;
       end

       if size(G.nodes.coords, 2) > 2,
          dz(k) = G.cells.volumes(ix(k))/(dx(k)*dy(k));
       else
          dz(k) = 0;
       end
    end
end




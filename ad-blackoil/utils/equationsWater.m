function [problem, state] = equationsWater(state0, state, model, dt, drivingForces, varargin)
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
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
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

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'iteration', -1);  % Compatibility only

opt = merge_options(opt, varargin{:});

W = drivingForces.W;

s = model.operators;
G = model.G;
f = model.fluid;

[p, wellSol] = model.getProps(state, 'pressure', 'wellsol');

[p0] = model.getProps(state0, 'pressure');

pBH    = vertcat(wellSol.bhp);
qWs    = vertcat(wellSol.qWs);

%Initialization of independent variables ----------------------------------

if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [p, qWs, pBH] = ...
            initVariablesADI(p, qWs, pBH);
    else
        [p0, tmp, tmp] = ...
            initVariablesADI(p0, sW0,          ...
            zeros(size(qWs)), ...
            zeros(size(qOs)), ...
            zeros(size(pBH)));                          %#ok
    end
end
primaryVars = {'pressure','qWs', 'bhp'};

clear tmp
%grav  = gravity;
%gdz   = s.Grad(G.cells.centroids) * model.gravity';
gdz   = s.Grad(G.cells.centroids) * model.getGravityVector()';
%--------------------
%check for p-dependent tran mult:
trMult = 1;
if isfield(f, 'tranMultR'), trMult = f.tranMultR(p); end

%check for p-dependent porv mult:
pvMult = 1; pvMult0 = 1;
if isfield(f, 'pvMultR')
    pvMult =  f.pvMultR(p);
    pvMult0 = f.pvMultR(p0);
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
upcw = (double(dpW)<=0);
vW = - s.faceUpstr(upcw, mobW).*trans.*dpW;
bWvW = s.faceUpstr(upcw, bW).*vW;

if model.outputFluxes
    state = model.storeFluxes(state, vW, [], []);
end

if model.extraStateOutput
    state = model.storebfactors(state, bW, [], []);
    state = model.storeMobilities(state, mobW, [], []);
    state = model.storeUpstreamIndices(state, upcw, [], []);
end
% EQUATIONS ---------------------------------------------------------------
% water:
eqs{1} = (s.pv/dt).*( pvMult.*bW - pvMult0.*f.bW(p0) ) + s.Div(bWvW);
eqs = addFluxesFromSourcesAndBC(model, eqs, ...
                                       {p},...
                                       {rhoW},...
                                       {mobW}, ...
                                       {bW},  ...
                                       {ones(numel(p0),1)}, ...
                                       drivingForces);

names = {'water'};
types = {'cell'};
% well equations
if ~isempty(W)
    if ~opt.reverseMode
        wc    = vertcat(W.cells);
        pw   = p(wc);
        rhos = [f.rhoWS];
        bw   = {bW(wc)};
        mw   = {mobW(wc)};
        s = {1};

        wm = WellModel();
        [cqs, weqs, ctrleqs, wc, state.wellSol]  = wm.computeWellFlux(model, W, wellSol, ...
                                             pBH, {qWs}, pw, rhos, bw, mw, s, {},...
                                             'nonlinearIteration', opt.iteration,'referencePressureIndex', 1);
        eqs(2) = weqs;
        eqs{3} = ctrleqs;
        
        eqs{1}(wc) = eqs{1}(wc) - cqs{1};        
        names(2:3) = {'waterWells', 'closureWells'};
        types(2:3) = {'perf', 'well'};
    else
        % in reverse mode just gather zero-eqs of correct size
        for eqn = 2:3
            nw = numel(state0.wellSol);
            zw = double2ADI(zeros(nw,1), p0);
            eqs(2:3) = {zw, zw};
        end
        names(2:3) = {'empty', 'empty'};
        types(2:3) = {'none', 'none'};
    end
end
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end
%--------------------------------------------------------------------------

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







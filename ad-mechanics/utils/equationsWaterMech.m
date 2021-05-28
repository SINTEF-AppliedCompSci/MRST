function [eqs, names, types, state] = equationsWaterMech(p0, state0, p, wellVars, ...
                                                      state, model, dt, mechTerm, ...
                                                      drivingForces, varargin)
%
%
% SYNOPSIS:
%   function [eqs, names, types, state] = equationsWaterMech(p0, state0, p, wellVars, state, model, dt, mechTerm, drivingForces, varargin)
%
% DESCRIPTION:
%
% PARAMETERS:
%   p0            - pressure (previous time step)
%   state0        - state    (previous time step)
%   p             - pressure
%   wellVars      - Well variables
%   state         - state (current time step)
%   model         - model class instance that is used.
%   dt            - time step
%   mechTerm      - mechanical input which will enter the computation of the
%                   effective porosity
%   drivingForces - structure that gathers the well parameters and boundary conditions.
%   varargin      - 
%
% RETURNS:
%   eqs   - The residual values as ADI variables (that is with the Jacobian)
%           if the inputs were also ADI.
%   names - The name of each equations
%   types - The type of each equations
%   state - Some field related to well control of the state variables may be updated.
%
% EXAMPLE:
%   run2DCase, runNorneExample
%
% SEE ALSO:
%   equationsWater

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

    % Note that state is given only for output
    opt = struct('iteration', -1, ...
                 'resOnly', false); % just to avoid warning
    opt = merge_options(opt, varargin{:});

    % Shorter names for some commonly used parts of the model and forces.
    s = model.operators;
    f = model.fluid;
    rock = model.rock;
    G = model.G;
    W = drivingForces.W;


    % grav  = gravity;
    % gdz   = s.Grad(G.cells.centroids) * model.gravity';
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
    bW   = f.bW(p);
    rhoW = bW.*f.rhoWS;
    % rhoW on face, avarge of neighboring cells (E100, not E300)
    rhoWf = s.faceAvg(rhoW);
    mobW  = trMult./f.muW(p);
    dpW   = s.Grad(p) - rhoWf.*gdz;
    % water upstream-index
    upcw = (value(dpW)<=0);
    vW   = - s.faceUpstr(upcw, mobW).*trans.*dpW;
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

    effPorVol = G.cells.volumes.*(rock.poro.*pvMult + rock.alpha .* mechTerm.new);
    effPorVol0 = G.cells.volumes.*(rock.poro.*pvMult0 + rock.alpha .* mechTerm.old);

    eqs{1} = (1 ./ dt) .* (effPorVol .* bW - effPorVol0.* f.bW(p0)) + s.Div(bWvW);

    [eqs, ~, qRes] = addFluxesFromSourcesAndBC(model, eqs, {p}, {rhoW}, {mobW}, {1}, drivingForces);

    if model.outputFluxes
       state = model.storeBoundaryFluxes(state, qRes, [], [], drivingForces);
    end
    
    names = {'water'};
    types = {'cell'};

    % Finally, add in and setup well equations
    wellSol = model.getProp(state, 'wellsol');
    [~, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);
    wellSol0 = model.getProps(state0, 'wellSol');

    [eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, ...
                                                      names, types, wellSol0, ...
                                                      wellSol, wellVars, ...
                                                      wellMap, p, {mobW}, ...
                                                      {rhoW}, {}, {}, dt, ...
                                                      opt);

end

function state = updateStateBlackOilGeneric(model, state, problem, dx, drivingForces)
%Generic update function for blackoil-like models
%
% SYNOPSIS:
%   state = updateStateBlackOilGeneric(model, state, problem, dx)
%
% DESCRIPTION:
%   This is a relatively generic update function that can dynamically work
%   out where increments should go based on the model implementation. It
%   can be used for simple models or used as inspiration for more exotic
%   models.
%
%   Presently handles either 2/3-phase with disgas/vapoil or n-phase
%   without dissolution.
%
% REQUIRED PARAMETERS:
%   model   - PhysicalModel subclass.
%
%   state   - State which is to be updated.
%
%   problem - Linearized problem from which increments were obtained
%
%   dx      - Increments created by solving the linearized problem.
%
%   drivingForces - Wells etc.
%
% RETURNS:
%   state - Updated state.

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

state0 = state;
W = drivingForces.Wells;

[disgas, vapoil] = deal(false);

if isprop(model, 'vapoil')
    vapoil = model.vapoil;
end

if isprop(model, 'disgas')
    disgas = model.disgas;
end
[restVars, satVars, wellVars] = model.splitPrimaryVariables(problem.primaryVariables); 


state = model.updateStateFromIncrement(state, dx, problem,...
                                            'pressure', model.dpMax);
restVars = model.stripVars(restVars, 'pressure');

saturations = lower(model.saturationVarNames);
% satSolVar = intersect(lower(problem.primaryVariables), saturations);


if (disgas || vapoil)
    % Black oil with dissolution
    so = model.getProp(state, 'so');
    sw = model.getProp(state, 'sw');
    sg = model.getProp(state, 'sg');
    
    % Magic status flag, see inside for doc
    st = getCellStatus(state0, so, sw, sg, disgas, vapoil);

    dr = model.getIncrement(dx, problem, 'x');
    dsw = model.getIncrement(dx, problem, 'sw');
    % Interpretation of "gas" phase varies from cell to cell, remove
    % everything that isn't sG updates
    dsg = st{3}.*dr - st{2}.*dsw;

    if disgas
        state = model.updateStateFromIncrement(state, st{1}.*dr, problem, 'rs', model.drsMax);
    end
    
    if vapoil
        state = model.updateStateFromIncrement(state, st{2}.*dr, problem, 'rv', model.drsMax);
    end

    dso = -(dsg + dsw);

    ds = zeros(numel(so), numel(saturations));
    ds(:, strcmpi(saturations, 'sw')) = dsw;
    ds(:, strcmpi(saturations, 'so')) = dso;
    ds(:, strcmpi(saturations, 'sg')) = dsg;
    
    state = model.updateStateFromIncrement(state, ds, problem, 's', inf, model.dsMax);
    % We should *NOT* be solving for oil saturation for this to make sense
    assert(~any(strcmpi(satVars, 'so')));
    state = computeFlashBlackOil(state, state0, model, st);
    state.s  = bsxfun(@rdivide, state.s, sum(state.s, 2));
    
    %  We have explicitly dealt with rs/rv properties, remove from list
    %  meant for autoupdate.
    restVars = model.stripVars(restVars, {'rs', 'rv', 'x'});
else
    state = model.updateSaturations(state, dx, problem, satVars);
end

% Update components
for i = 1:numel(restVars);
     p = restVars{i};
     % Update the state
     state = model.updateStateFromIncrement(state, dx, problem, p);
end

% Update the wells
state.wellSol = model.updateWellSol(state.wellSol, dx, problem, wellVars);


% Handle the directly assigned values (i.e. can be deduced directly from
% the well controls.
wi = strcmpi(saturations, 'sw');
oi = strcmpi(saturations, 'so');
gi = strcmpi(saturations, 'sg');

state.wellSol = assignWellValuesFromControl(model, state.wellSol, W, wi, oi, gi);
end

function st = getCellStatus(state, oil, wat, gas, disgas, vapoil)
% Status should be passed on from updateStateVO (to be sure definition is
% identical). rs and rv are assumed to be compatible, i.e. rx = rxSat for
% saturated cells and rx <= rxSat for undersaturated. Three values of
% status are:
% status 0: should not occur (almost water only -> state 3)
% status 1 oil, no gas  : x = rs, sg = 0    , rv = rvMax
% status 2 gas, no oil  : x = rv, sg = 1-sw , rs = rsMax
% status 3 oil and gas  : x = sg, rs = rsMax, rv = rvMax
if isfield(state, 'status')
    status = state.status;
else
    watOnly    = wat > 1- sqrt(eps);
    if ~vapoil
        oilPresent = true;
    else
        oilPresent = or(oil > 0, watOnly);
    end
    if ~disgas
        gasPresent = true;
    else
        gasPresent = or(gas > 0, watOnly);
    end
    status = oilPresent + 2*gasPresent;
end

if ~disgas
    st1 = false;
else
    st1 = status==1;
end
if ~vapoil
    st2 = false;
else
    st2 = status==2;
end
st3 = status == 3;
st = {st1, st2, st3};
end


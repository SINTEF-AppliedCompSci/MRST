function state = computeFlashBlackOil(state, state0, model, status)
% Compute flash for a black-oil model with disgas/vapoil
%
% SYNOPSIS:
%   state = computeFlashBlackOil(state, state0, model, status)
%
% DESCRIPTION:
%   Compute flash to ensure that dissolved properties are within physically
%   reasonable values, and simultanously avoid that properties go far
%   beyond the saturated zone when they were initially unsaturated and vice
%   versa.
%
% REQUIRED PARAMETERS:
%   state  - State where saturations, rs, rv have been updated due to a
%            Newton-step.
%
%   state0 - State from before the linearized update.
%
%   model  - The ThreePhaseBlackOil derived model used to compute the
%            update.
%
%   status - Status flags from getCellStatusVO, applied to state0.
%
% RETURNS:
%   state - Updated state where saturations and values are chopped near
%           phase transitions.
%
% NOTE:
%    Be mindful of the definition of state0. It is not necessarily the
%    state at the previous timestep, but rather the state at the previous
%    nonlinear iteration!

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
    disgas = model.disgas;
    vapoil = model.vapoil;

    so = model.getProp(state, 'so');
    so0 = model.getProp(state0, 'so');
    
    if model.water
        sw = model.getProp(state, 'sw');
    else
        sw = zeros(model.G.cells.num, 1);
    end
    
    sg = model.getProp(state, 'sg');
    sg0 = model.getProp(state0, 'sg');
    
    rs = model.getProp(state, 'rs');
    rs0 = model.getProp(state0, 'rs');
    
    rv = model.getProp(state, 'rv');
    rv0 = model.getProp(state0, 'rv');
    
    etol = sqrt(eps);
    % Determine status of updated cells -----------------------------------------
    watOnly  = sw > 1-etol;
    % phase transitions sg <-> rs  --------------------------------------------
    if ~disgas
        rsSat0 = rs0;
        rsSat  = rsSat0; 
        gasPresent = true;
    else
        st1 = status{1};
        pvt = model.PVTPropertyFunctions;
        rsSat0 = pvt.get(model, state0, 'RsMax');
        rsSat = pvt.get(model, state, 'RsMax');
        gasPresent = or(and( sg > 0 | rs == 0, ~st1), watOnly); % Obvious case
        % Keep oil saturated if previous sg is sufficiently large:
        ix1 = and( sg < 0, sg0 > etol);
        gasPresent = or(gasPresent, ix1);
        % Set oil saturated if previous rs is sufficiently large
        ix2 = and( and(rs > rsSat*(1+etol), st1), rs0 > rsSat0*(1-etol) );
        sg(ix2) = 0;
        gasPresent = or(gasPresent, ix2);
    end
    ix = sg < 0;
    sw(ix) = sw(ix)./(1-sg(ix));
    so(ix) = so(ix)./(1-sg(ix));
    sg(ix) = 0;

    % phase transitions so <-> rv
    if ~vapoil
        oilPresent = true;
        rvSat0 = rv0;
        rvSat  = rvSat0;
    else
        st2 = status{2};
        rvSat0 = model.getProp(state0, 'RvMax');
        rvSat = model.getProp(state, 'RvMax');
        oilPresent = or(and( so > 0 | rv == 0, ~st2), watOnly); % Obvious case
        % Keep gas saturated if previous so is sufficiently large
        ix1 = and( so < 0, so0 > etol);
        oilPresent = or(oilPresent, ix1);
        % Set gas saturated if previous rv is sufficiently large
        ix2 = and( and(rv > rvSat*(1+etol), st2), rv0 > rvSat0*(1-etol) );
        so(ix2) = 0;
        oilPresent = or(oilPresent, ix2);
    end
    ix = so < 0;
    sw(ix) = sw(ix)./(1-so(ix));
    sg(ix) = sg(ix)./(1-so(ix));
    so(ix) = 0;

    % make sure sw >=0
    ix = sw < 0;
    so(ix) = so(ix)./(1-sw(ix));
    sg(ix) = sg(ix)./(1-sw(ix));
    sw(ix) = 0;

    % Update saturated r-values -----------------------------------------------
    rs(gasPresent) = rsSat(gasPresent);
    rv(oilPresent) = rvSat(oilPresent);

    % Update undersatured r-values
    rs(~gasPresent) = min(rsSat(~gasPresent), rs(~gasPresent));

    % Update state ------------------------------------------------------------
    if model.water
        state = model.setProp(state, 'sw', sw);
    end
    state = model.setProp(state, 'so', so);
    state = model.setProp(state, 'sg', sg);
    
    state = model.setProp(state, 'rs', max(rs, 0));
    state = model.setProp(state, 'rv', max(rv, 0));

    state.status = oilPresent + 2*gasPresent;
end


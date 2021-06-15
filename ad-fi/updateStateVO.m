function [state] = updateStateVO(W, state, dx, f, system)
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

stepOpts = system.stepOptions;

disgas = system.activeComponents.disgas;
vapoil = system.activeComponents.vapoil;
[st1, st2, st3] = getCellStatus(state, disgas, vapoil);
etol = sqrt(eps);

dpMax = stepOpts.dpMax;
dp    = dx{1};
dp    = sign(dp).*min(abs(dp), abs(dpMax.*state.pressure));
p     = state.pressure + dp;
p     = max(1*barsa, p);

% saturation updates
dsw = dx{2};
dsg = st3.*dx{3} - st2.*dsw;
dso = -dsw-dsg;
maxVal = max(abs([dsw, dso, dsg]), [], 2);
step   = min(stepOpts.dsMax./maxVal, 1);
sw = state.s(:,1) + step.*dsw;
so = state.s(:,2) + step.*dso;
sg = state.s(:,3) + step.*dsg;

% r updates
drMax = stepOpts.drsMax;
if disgas
    rs  = state.rs;
    drs = st1.*dx{3};
    drs = sign(drs).*min(abs(drs), abs(drMax.*rs));
    rs  = rs + drs;
end
if vapoil
    rv  = state.rv;
    drv = st2.*dx{3};
    drv = sign(drv).*min(abs(drv), abs(drMax.*rv));
    rv  = rv + drv;
end

% Determine status of updated cells -----------------------------------------
watOnly  = sw(:,1) > 1-etol;

% phase transitions sg <-> rs  --------------------------------------------
if ~disgas
    rsSat0 = state.rs;
    rsSat  = rsSat0; 
    gasPresent = true;
else
    rsSat0 = f.rsSat(state.pressure);
    rsSat  = f.rsSat(p);
    gasPresent = or(and( sg > 0, ~st1), watOnly); % Obvious case
    % Keep oil saturated if previous sg is sufficiently large:
    ix1 = and( sg < 0, state.s(:,3) > etol);
    gasPresent = or(gasPresent, ix1);
    % Set oil saturated if previous rs is sufficiently large
    ix2 = and( and(rs > rsSat*(1+etol), st1), state.rs > rsSat0*(1-etol) );
    assert(all(sg(ix2)==0))
    gasPresent = or(gasPresent, ix2);
end
ix = sg < 0;
sw(ix) = sw(ix)./(1-sg(ix));
so(ix) = so(ix)./(1-sg(ix));
sg(ix) = 0;

% phase transitions so <-> rv
if ~vapoil
    oilPresent = true;
    rvSat0 = state.rv;
    rvSat  = rvSat0;
else
    rvSat0   = f.rvSat(state.pressure);
    rvSat    = f.rvSat(p);
    oilPresent = or(and( so > 0, ~st2), watOnly); % Obvious case
    % Keep gas saturated if previous so is sufficiently large
    ix1 = and( so < 0, state.s(:,2) > etol);
    oilPresent = or(oilPresent, ix1);
    % Set gas saturated if previous rv is sufficiently large
    ix2 = and( and(rv > rvSat*(1+etol), st2), state.rv > rvSat0*(1-etol) );
    assert(all(so(ix2)==0))
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
state.pressure = p;
% saturations should now be ok, but make proj to be sure
state.s  = bsxfun(@rdivide, [sw so sg], sw + so + sg);
state.rs = max(rs, 0);
state.rv = max(rv, 0);
state.status = oilPresent + 2*gasPresent;

% Wells -------------------------------------------------------------------
dqWs = dx{4};
dqOs = dx{5};
dqGs = dx{6};
dpBH = dx{7};

if ~isempty(dpBH)  % Avoid case where there is no active well.
   dpBH = sign(dpBH).*min(abs(dpBH), abs(dpMax.*vertcat(state.wellSol.bhp)));
   for w = 1:numel(state.wellSol)
      ws = state.wellSol(w);
      ws.bhp  = ws.bhp + dpBH(w);
      ws.qWs  = ws.qWs + dqWs(w);
      ws.qOs  = ws.qOs + dqOs(w);
      ws.qGs  = ws.qGs + dqGs(w);

      tp = ws.type;
      v  = ws.val;
      switch tp
        case 'bhp'
          ws.bhp = v;
        case 'rate'
          ws.qWs = v*W(w).compi(1);
          ws.qOs = v*W(w).compi(2);
          ws.qGs = v*W(w).compi(3);
        case 'orat'
          ws.qOs = v;
        case 'wrat'
          ws.qWs = v;
        case 'grat'
          ws.qGs = v;
      end
      state.wellSol(w) = ws;
   end
end

end

function [st1, st2, st3] = getCellStatus(state, disgas, vapoil)
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
    s = state.s;
    watOnly    = s(:,1) > 1- sqrt(eps);
    if ~vapoil
        oilPresent = true;
    else
        oilPresent = or(s(:,2) > 0, watOnly);
    end
    if ~disgas
        gasPresent = true;
    else
        gasPresent = or(s(:,3) > 0, watOnly);
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
end

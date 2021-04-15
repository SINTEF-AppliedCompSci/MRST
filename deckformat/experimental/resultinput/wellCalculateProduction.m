function report = wellCalculateProduction(state, W, fluid, time)
%Convert MRST well solution data to ECLIPSE keyword representation.
%
% SYNOPSIS:
%   report = wellCalculateProduction(rSol, wSol, W, fluid, time)
%
% PARAMETERS:
%   rSol  - Reservoir solution structure as defined by a flow and/or a
%           transport solver routine
%
%   wSol  - Well solution data structure.
%
%   W     - Well data structure as defined by function 'addWell'.
%           Report data will be generated for each well in 'W'.
%
%   fluid - Fluid data structure.
%
%   time  - Report time.
%
% RETURNS:
%   report - ECLIPSE keyword representation of the well solution data.
%            A structure array containing the following fields:
%              - TIME -- Report time (== time).
%              - WBHP -- Bottom-hole pressure in each well at TIME.
%              - WVPT -- Total rate (sum of all perforation fluxes).
%              - WOPR -- Total oil production.
%              - WWPR -- Total water production.
%              - WWCT -- Well water-cut.
%
% NOTE:
%   Multiple 'report' structures--e.g., from several time steps--may be
%   concatenated using the 'addToTimeStruct' function.
%
% SEE ALSO:
%   `incompTPFA`, `explicitTransport`, `addWell`, `addToTimeStruct`.

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


   if isempty(W),
      report = struct('TIME', [], 'WBHP', [], 'WVPT', [], ...
                      'WOPR', [], 'WWPR', [], 'WWCT', []);
   else
      arg = { 'UniformOutput', false };
      wSol = state.wellSol;

      mu = fluid.properties(state);
      s  = fluid.saturation(state);
      kr = fluid.relperm(s(:,1), state);

      mob = bsxfun(@rdivide, kr, mu);
      f   = bsxfun(@rdivide, mob, sum(mob, 2));
      fw  = cellfun(@(c) f(c,1), { W.cells }, arg{:});

      rate0 = @(flx, f) sum([flx, flx.*(1-f), flx.*f], 1);
      rate1 = @(f,w) rate0(reshape(f, [], 1), reshape(w, [], 1));
      rates = cellfun(rate1, { wSol.flux }, fw, 'UniformOutput', false);
      rates = vertcat(rates{:});  % All 'rates' elements are 1-by-3 DOUBLE.

      WBHP = [ wSol.pressure ] .';
      WVPT = rates(:,1);
      WOPR = rates(:,2);
      WWPR = rates(:,3);
      WWCT = WWPR ./ WVPT;

      report = struct('TIME', time, 'WBHP', WBHP, 'WVPT', WVPT, ...
                      'WOPR', WOPR, 'WWPR', WWPR, 'WWCT', WWCT);
   end
end

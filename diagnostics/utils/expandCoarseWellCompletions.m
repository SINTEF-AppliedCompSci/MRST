function [exc, eWC] = expandCoarseWellCompletions(xc, WC, Wdf, p)
%Pseudo-wells for computing flow diagnostics in an upscaled model
%
% SYNOPSIS:
%   [xn, Wn] = expandCoarseWellCompletions(xc, Wc, Wdf, p)
%
% REQUIRED PARAMETERS:
%   xc  - Reservoir and well solution structure either properly
%         initialized from functions 'initResSol' and 'initWellSol',
%         respectively, or the results from a call to function
%         'solveIncompFlow'.  Must contain valid cell interface fluxes,
%         'state.flux'.
%
%   Wc  - Well structure as defined by function 'addWell'.
%
%   Wdf - Expanded well structure used for flow diagnostics computation in
%         the underlying fine-scale model, typically the result of a call
%         to 'expandWellCompletions'.
%
%   p   - partition vector describing the coarse-scale partition (i.e., the
%         mapping from the fine grid to the coarse grid)
%
% RETURNS:
%   Wn  - Well structure containing the pseudo wells
%
%   xn  - Reservoir and solution structure corresponding
%
% SEE ALSO:
%   `expandWellCompletions`, `computeTOFandTracer`, `computeWellPairs`

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

eWC = [];
exc = xc;
mp  = max(p);
n=0;
for i=1:numel(Wdf)
   for j=1:numel(WC)
      flag = false(mp,1);
      flag(unique(p(vertcat(Wdf(i).cells))))=true;
      ind = flag(WC(j).cells);
      if ~sum(ind), continue, end;
      
      nW = WC(j);
      nW.name = Wdf(i).name;
      nW.cells = WC(j).cells(ind);
      if numel(WC(j).r)==1
         nW.r = WC(j).r;
      else
         nW.r = WC(j).r(ind);
      end
      nW.dir = WC(j).dir(ind);
      nW.WI  = WC(j).WI(ind);
      nW.dZ  = WC(j).dZ(ind);
      eWC = [eWC; nW];                                          %#ok<AGROW>
      n = n+1;
      exc.wellSol(n).flux = xc.wellSol(j).flux(ind);
   end
end
end


function [vpt, opr, wpr, wc] = prodCurves(w, state, fluid)
%Calculate derived production data for wells.
%
% SYNOPSIS:
%   [tot, opr, wpr, wc] = prodCurves(W, state, fluid)
%
% PARAMETERS:
%   W     - Well data structure as defined by function addWell.
%
%   state - Solution state object as defined by functions 'initResSol',
%           'initWellSol' and, possibly, flow and/or transport solvers.
%
%   fluid - Fluid data structure.
%
% RETURNS:
%   tot - Total well production.  One value for each individual well.
%
%   opr - Oil (liquid) production rate.  One value for each individual well.
%
%   wpr - Water production rate.  One value for each individual well.
%
%   wc  - Water cut.  One value for each individual well.
%
% LIMITATION:
%   This function assumes that the solution state object represents
%   incompressible two-phase flow.

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


   arg = { 'UniformOutput', false };

   mu = fluid.properties(state);
   s  = fluid.saturation(state);
   kr = fluid.relperm(s(:,1), state);

   mob = bsxfun(@rdivide, kr, mu);
   f   = bsxfun(@rdivide, mob, sum(mob, 2));
   fw  = cellfun(@(c) f(c,1), { w.cells }, arg{:});

   rate0 = @(flx, f) sum([flx, flx.*(1-f), flx.*f], 1);
   rate1 = @(flx, f) rate0(reshape(-flx, [], 1), reshape(f, [], 1));
   rates = cellfun(rate1, { state.wellSol.flux }, fw, arg{:});
   rates = vertcat(rates{:});

   vpt = rates(:,1);
   opr = rates(:,2);
   wpr = rates(:,3);
   wc  = wpr ./ vpt;
end

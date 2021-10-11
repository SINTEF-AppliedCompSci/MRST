function dp = computeIncompWellPressureDrop(W, mob, rho, g)
%Compute incompressible connection pressure drop for a single well
%
% SYNOPSIS:
%   dp = computeIncompWellPressureDrop(W, mob, rho, g)
%
% PARAMETERS:
%   W   - Well structure as defined by function addWell.
%
%   mob - Phase mobilities for each cell in model.
%
%   rho - Fluid mass densities.  One scalar for each fluid phase.
%
%   g   - Norm of gravity along z-axis
%
% RETURNS:
%   dp - Column vector of size numel(W.cells).  Contains each completion's
%        pressure drop along the well bore from the bottom hole pressure
%        reference depth, under certain assumptions for the incompressible,
%        linear pressure equation.
%
% NOTES:
%   In order to avoid nonlinear behavior for wells, this function assumes
%   that the well bore density is constant for all perforations and is
%   equal to the injection composition if the well is an injector,
%   either through positive rates or W.sign > 0. If the well is a producer,
%   the mixture is assumed to be taken from the perforated cells, with flux
%   proportional to the well indices.
%
% SEE ALSO:
%   `addWell`, `incompTPFA`, `incompMimetic`, `incompMPFA`

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

   check_input(W, mob, rho)

   np = numel(rho);

   if is_injector(W)
      % We compute a simple mixture of the injected fluids
      i = 1 : min(np, size(W.compi, 2));

      rhoMix = rho(i) * W.compi(i)';
   else
      % We compute approximate phase fluxes by taking the fractional flow
      % in each cell, weighted by the densities.  The flux is estimated
      % roughly as proportional to the well index.
      i = 1 : min(np, size(mob, 2));

      mobw = mob(W.cells, :);
      f    = bsxfun(@rdivide, mobw, sum(mobw, 2));

      rhoMix = sum(W.WI .* (f(:, i) * rho(i)')) ./ sum(W.WI);
   end

   dp = g * (W.dZ * rhoMix);
end

%--------------------------------------------------------------------------

function check_input(W, mob, rho)
   persistent NI NP NM

   nchanged = 0;

   if isempty(NI), NI = numel(W.compi); nchanged = nchanged + 1; end
   if isempty(NP), NP = numel(rho);     nchanged = nchanged + 1; end
   if isempty(NM), NM = size(mob, 2);   nchanged = nchanged + 1; end

   [NP, nchanged] = change_if_modified(NP, nchanged, numel(rho));
   [NI, nchanged] = change_if_modified(NI, nchanged, numel(W.compi));
   [NM, nchanged] = change_if_modified(NM, nchanged, size(mob, 2));

   if (nchanged > 0) && ((NI ~= NP) || (NM ~= NP))
      np = min([NI, NM, NP]);
      pl = '';  if np ~= 1, pl = 's'; end

      warning('NumPhase:Inconsistent', ...
             ['Inconsistent Number of Phases.  ', ...
              'Using %d Phase%s (=min([%d, %d, %d])).'], ...
              np, pl, NI, NP, NM);
   end
end

%--------------------------------------------------------------------------

function tf = is_injector(W)
   % A well is an injector if the sign is set (and positive) or the well is
   % constrained to a positive rate target.  Otherwise, we assume that the
   % well is a producer.
   tf = (W.sign > 0) || (strcmpi(W.type, 'rate') && (W.val > 0));
end

%--------------------------------------------------------------------------

function [n, count] = change_if_modified(n, count, new)
   if n ~= new
      n     = new;
      count = count + 1;
   end
end

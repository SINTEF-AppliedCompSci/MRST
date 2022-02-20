function WP = computeWellPairs(state, G, rock, W, D )
%Compute volumes and fluxes associated with each flux pair
%
% SYNOPSIS:
%   WP = computeWellPairs(state, G, rock, W, D)
%
% PARAMETERS:
%   state - Reservoir state
%
%   G     - Grid structure.
%
%   rock  - Rock data structure.
%           Must contain a valid porosity field, 'rock.poro'.
%
%   state - Reservoir and well solution structure either properly
%           initialized from functions 'initResSol' and 'initWellSol'
%           respectively, or the results from a call to function
%           'solveIncompFlow'.  Must contain valid cell interface fluxes,
%           'state.flux'.
%
%   W     - Well structure
%
%   D     - Struct containing TOF and tracer
%           See documentation of computeTOFandTracer
% RETURNS:
%   WP - struct containing well-pair diagnostics
%       'pairs'   - list of well pairs (as text strings using well names)
%       'pairsIx' - list of well pairs (indices in a N_pair x 2 list. Pair
%                   number i has injector pairsIx(i, 1) and producer
%                   pairsIx(i, 2)
%       'vols'    - pore volumes associated with each pair
%       'inj'     - list of structs containing allocation factors for the
%                   injectors
%       'prod'    - list of structs containing allocation factors for the
%                   producers
%
%   The structs that contain flux allocation information consist of the
%   following data objects:
%       'alloc'   - nxm array, where n is the number of perforations in
%                   this well and m is the number of wells or segments from
%                   the flow diagnostics computation that may potentially
%                   contribute to the flux allocation
%       'ralloc'  - nx1 array containing the well flux that cannot be
%                   attributed to one of the wells/segments accounted for
%                   in the flow diagnostics computatation
%       'z'       - nx1 array giving the depth of the completions
%       'name'    - string with the name of the well or segment

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


if ~isfield(state.wellSol, 'flux') && isfield(state.wellSol, 'cqs')
   for i=1:numel(state.wellSol)
      state.wellSol(i).flux = sum(state.wellSol(i).cqs, 2);
   end
end

ni = numel(D.inj);
np = numel(D.prod);

% Names of well pairs
WP.pairs = cellfun(@(i,p) [i, ', ', p], ...
   repmat({W(D.inj).name},[1,np]), rldecode({W(D.prod).name},ni,2), ...
   'UniformOutput', false);

vols = D.itracer'*bsxfun(@times, poreVolume(G, rock), D.ptracer);
% Volumes associated with each well pair
WP.vols = vols(:)';

% Numerical indices of each pair
WP.pairIx = [repmat(1:ni, 1, np);...
             reshape(repmat(1:np, ni, 1), 1, [])]';
% Compute allocation factors
for i=1:ni
   qik = sum(state.wellSol(D.inj(i)).flux, 2);
   if isempty(qik)
       continue;
   end
   cj  = D.ptracer(W(D.inj(i)).cells,:);
   WP.inj(i).alloc  = repmat(qik, 1, np) .* cj;% / sum(qik);
   WP.inj(i).ralloc = qik - sum(WP.inj(i).alloc, 2);
   WP.inj(i).z      = W(D.inj(i)).dZ + W(D.inj(i)).refDepth;
   WP.inj(i).name   = W(D.inj(i)).name;
end 
for i=1:np
   qjk = sum(state.wellSol(D.prod(i)).flux, 2);
   if isempty(qjk)
       continue;
   end
   ci  = D.itracer(W(D.prod(i)).cells,:);
   WP.prod(i).alloc  = repmat(qjk, 1, ni) .* ci;% / sum(qjk);
   WP.prod(i).ralloc = qjk - sum(WP.prod(i).alloc, 2);
   WP.prod(i).z      = W(D.prod(i)).dZ + W(D.prod(i)).refDepth;
   WP.prod(i).name   = W(D.prod(i)).name;
end
end


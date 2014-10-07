function WP = computeWellPairs(state, G, rock, W, D )
%Compute volumes and fluxes associated with each flux pair
%
% SYNOPSIS:
%   WP = computeWellPairs(G, rock, W, D)
%
% PARAMETERS:
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
%       'pairs'   - list of well pairs
%       'vols'    - pore volumes associated with each pair

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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

% If only one injector / producer is present, ensure one over the domain
if ni == 1
    D.itracer = ones(G.cells.num, 1);
end

if np == 1
    D.ptracer = ones(G.cells.num, 1);
end

% Fraction of cell associated with each well pair
cij = repmat(D.itracer, 1, np) .* rldecode(D.ptracer, ni,2);

% Volumes associated with each well pair
WP.vols = sum(cij .* repmat(poreVolume(G, rock),1, ni*np));

% Compute allocation factors
for i=1:ni
   qik = state.wellSol(D.inj(i)).flux;
   ci  = D.itracer(W(D.inj(i)).cells,i);
   cj  = D.ptracer(W(D.inj(i)).cells,:);
   WP.inj(i).alloc = repmat(qik.*ci, 1, np) .* cj;% / sum(qik);
   WP.inj(i).z     = W(D.inj(i)).dZ;
end
for i=1:np
   qjk = state.wellSol(D.prod(i)).flux;
   ci  = D.itracer(W(D.prod(i)).cells,:);
   cj  = D.ptracer(W(D.prod(i)).cells,i);
   WP.prod(i).alloc = repmat(qjk.*cj, 1, ni) .* ci;% / sum(qjk);
   WP.prod(i).z     = W(D.prod(i)).dZ;
end
end


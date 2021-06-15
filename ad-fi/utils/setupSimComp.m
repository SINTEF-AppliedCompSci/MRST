function s = setupSimComp(G, rock, varargin)
% Set up helper structure for solvers based on automatic differentiation.
%
%
% SYNOPSIS:
%   s = setupSimComp(G, rock)
%
% PARAMETERS:
%   G       - Grid structure
%   rock    - rock structure
%
% OPTIONAL PARAMETERS
%   deck      - deck file containing rock properties
%   trans     - transmissibility for internal faces (if neighbors given) or for all faces (if
%               neighbors are not given)
%   neighbors - neighbors for each internal face
%   porv      - pore volumes for all cells

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

opt = struct('deck', [], 'neighbors', [], 'trans', [], 'porv', []);
opt = merge_options(opt, varargin{:});


T = opt.trans;
N = opt.neighbors;

if isempty(N)
   % Get neighbors for internal faces from grid.
   N  = double(G.faces.neighbors);
   intInx = (prod(N,2)~=0);
   N  = N(intInx, :);
else
   % neighbors are given for internal faces.
   intInx = 1:size(N, 1);
end

if isempty(T)
   % half-trans -> trans and reduce to interior
   m = [];
   if ~isempty(opt.deck)
      m = computeTranMult(G, opt.deck.GRID);
   end
   if isempty(m)
      m = 1;
   end
   T = m.*computeTrans(G, rock);
   cf = G.cells.faces(:,1);
   nf = G.faces.num;
   T  = 1 ./ accumarray(cf, 1./T, [nf, 1]);
   s.T_all = T;
   T = T(intInx);
else
   s.T_all = T;
   T = T(intInx);
end

s.T = T;

pv = opt.porv;
if isempty(pv)
    if isfield(G, 'PORV')
        pv = G.PORV;
    else
        pv = poreVolume(G, rock);
    end
end
s.pv = pv;

% C - (transpose) divergence matrix
n = size(N,1);
C  = sparse( [(1:n)'; (1:n)'], N, ones(n,1)*[1 -1], n, G.cells.num);
s.C = C;
s.grad = @(x)C*x;
s.div  = @(x)C'*x;

% faceAvg - as multiplication with matrix
nc = max(max(N));
nf = size(N,1);
M  = sparse((1:nf)'*[1 1], N, .5*ones(nf,2), nf, nc);
s.faceAvg = @(x)M*x;

% faceUpstr - as multiplication with matrix
s.faceUpstr = @(flag, x)faceUpstr(flag, x, N, [nf, nc]);

% Include neighbor relations since eqsfiOW needs it...
s.N = N;
end

function xu = faceUpstr(flag, x, N, sz)
upCell       = N(:,2);
upCell(flag) = N(flag,1);
xu = sparse((1:sz(1))', upCell, 1, sz(1), sz(2))*x;
end



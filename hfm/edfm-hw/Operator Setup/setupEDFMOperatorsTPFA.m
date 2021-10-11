function s = setupEDFMOperatorsTPFA(G, rock, tol, varargin)
% Modified from setupOperatorsTPFA.m. Additional lines are added to include
% non-neighbouring connections. Do not use the optional parameters when
% using this function. See Lines 109 and 133.
%
% Set up helper structure for solvers based on automatic differentiation.
%
% SYNOPSIS:
%   s = setupOperatorsTPFA(G, rock)
%
% DESCRIPTION:
%   The automatic differention solvers rely on discrete operators for
%   divergence and gradient on the grid as well as a variety of derived
%   reservoir quantities such as transmissibility and pore volume. The
%   purpose of this function is to assemble all such quantities using a
%   standard two-point finite volume approxiation (TPFA).
% 
% PARAMETERS:
%
%   G       - MRST grid given as a struct. See grid_structure.m for more
%             details.
%
%   rock    - Rock structure containing fields .perm and .poro with
%             approprioate dimensions for the grid. See makeRock for more
%             details.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%
%   'deck'      - deck file containing rock properties
%
%   'trans'     - transmissibility for internal faces (if neighbors given) or for all faces (if
%                 neighbors are not given)
%
%   'neighbors' - neighbors for each internal face
%
%   'porv'      - pore volumes for all cells
% RETURNS:
%   s           - Operators struct, with discrete operators and derived
%                 quantities:
%         T_all:  Transmissibilities for all interfaces, *including*
%                 (half) transmissibilities for faces on the boundary. One
%                 value per interface.
%
%         T    : Transmissibilities for all internal interfaces. Internal
%                interfaces have a cell on both sides.
%
%         pv   : Pore volumes. See function poreVolume. One value per cell.
%
%         C    : Transfer matrix between cells and faces. Used to derive
%                discrete gradient and divergence operators.
%
%      faceAvg : (Function) For each interface, computes the average value
%      of a quantity defined in the cells. If a face is connecting cells 4
%      and 5, the faceAvg function will compute the arithmetic average of
%      the values in both cells.
%
%      Grad    : (Function)  Discrete gradient. Computes the gradient on
%      each interface via a first order finite difference approximation
%      using the values of the cells connected to the face.  
%                
%      Div     : (Function) Discrete divergence. Integrates / sums up
%      values on the interfaces for all cells to give the (integrated)
%      divergence per cell.
%
%      faceUpstr: (Function) Perform upstream weighting of values. Given a
%      set of cell wise values and a upstream flag for each interface, this
%      function will pick the values corresponding to the position in the
%      neighborship. I.e. if the flag is true for a given interface, it
%      will select the value in the FIRST cell connected to the interface
%      x(N(faceNo, 1)). Otherwise, it will select the SECOND x(N(faceNo,
%      2)). Typical usage is for upstream weighting of transported
%      quantities.
%
%      N      : Neighborship structure. Will be number of interfaces by 2
%               in size where N(ix, :) contains the cells connected to
%               face number ix.
%
%

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
   N=[N;G.nnc.cells]; % modified line!!!!!!!!!!!!
   intInx = all(N ~= 0, 2);
   N  = N(intInx, :);
else
   % neighbors are given
   n_if = size(N, 1);
   intInx = true(n_if, 1);
   if isfield(G, 'faces')
       % Try to match given interfaces to actual grid.
       intInxGrid = all(G.faces.neighbors ~= 0, 2);
       if sum(intInxGrid) == n_if
           % Given neighbors correspond to internal interfaces
           intInx = intInxGrid;
       elseif n_if == G.faces.num
           % Given neighbors correspond to *all* interfaces
           intInx = all(N ~= 0, 2);
       end
   end
end

if isempty(T)
   % half-trans -> trans and reduce to interior
   T = getFaceTransmissibility(G, rock, opt.deck);
   transmult=transmultEDFM(G,tol);
   T=T.*transmult; % multiply with transmissibility multiplier
   s.T_all = T; % T_all will not contain nnc transmissibilities!!!!!
   T=[T;G.nnc.T]; % modified line!!!!!!!!!!!!!!!!
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
    zeropv = find(pv == 0);
    if ~isempty(zeropv)
        warning(['I computed zero pore volumes in ', num2str(numel(zeropv)), ...
                 ' cells. Consider adjusting poro / ntg fields or grid.']);
    end
end
s.pv = pv;

% C - (transpose) divergence matrix
n = size(N,1);
C  = sparse( [(1:n)'; (1:n)'], N, ones(n,1)*[1 -1], n, G.cells.num);
s.C = C;
s.Grad = @(x) -C*x;
s.Div  = @(x) C'*x;

% faceAvg - as multiplication with matrix
nc = max(max(N));
nf = size(N,1);
M  = sparse((1:nf)'*[1 1], N, .5*ones(nf,2), nf, nc);
s.faceAvg = @(x)M*x;

% faceUpstr - as multiplication with matrix
upw = @(flag, x)faceUpstr(flag, x, N, [nf, nc]);
s.faceUpstr = upw;

s.splitFaceCellValue = @(operators, flag, x) splitFaceCellValue(operators, flag, x, [nf, nc]);

% Include neighbor relations
s.N = N;
s.internalConn = intInx;
end

function xu = faceUpstr(flag, x, N, sz)
    if numel(flag) == 1
        flag = repmat(flag, size(N, 1), 1);
    end
    assert(numel(flag) == size(N, 1) && islogical(flag), ...
        'One logical upstream flag must'' be supplied per interface.');
    upCell       = N(:,2);
    upCell(flag) = N(flag,1);
    if isnumeric(x)
        % x is a simple matrix, we just extract values using the cells
        xu = x(upCell, :);
    else
        % x is likely AD, construct a matrix to achieve upstream weighting
        xu = sparse((1:sz(1))', upCell, 1, sz(1), sz(2))*x;
    end
end



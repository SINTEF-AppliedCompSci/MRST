function s = setupOperatorsTPFA(G, rock, varargin)
% Set up helper structure for solvers based on automatic differentiation.
%
%
% SYNOPSIS:
%   s = setupOperatorsTPFA(G, rock)
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

opt = struct('deck', [], 'neighbors', [], 'trans', [], 'porv', []);
opt = merge_options(opt, varargin{:});


T = opt.trans;
N = opt.neighbors;

if isempty(N)
   % Get neighbors for internal faces from grid.
   N  = double(G.faces.neighbors);
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
s.Grad = @(x) -C*x;
s.Div  = @(x) C'*x;

% faceAvg - as multiplication with matrix
nc = max(max(N));
nf = size(N,1);
M  = sparse((1:nf)'*[1 1], N, .5*ones(nf,2), nf, nc);
s.faceAvg = @(x)M*x;

% faceUpstr - as multiplication with matrix
s.faceUpstr = @(flag, x)faceUpstr(flag, x, N, [nf, nc]);

% Include neighbor relations
s.N = N;
s.internalConn = intInx;

% Average face values to cells
fvg = 0.5*bsxfun(@rdivide, ones(n,2), G.faces.areas(intInx));
aC  = sparse( [(1:n)'; (1:n)'], N, fvg, n, G.cells.num);
s.cellAvg = @(x) aC*x;
end

function xu = faceUpstr(flag, x, N, sz)
    if numel(flag) == 1
        flag = repmat(flag, size(N, 1), 1);
    end
    assert(numel(flag) == size(N, 1) && islogical(flag), ...
        'One logical upstream flag must'' be supplied per interface.');
    upCell       = N(:,2);
    upCell(flag) = N(flag,1);
    xu = sparse((1:sz(1))', upCell, 1, sz(1), sz(2))*x;
end



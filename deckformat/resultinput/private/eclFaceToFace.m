function ix = eclFaceToFace(G, NNC, varargin)
% ix = eclFaceToFace(G, NNC)
% Produce indices of interior faces such that an eclipse face propery F 
% converts to grid-face property f by
%  f(iI) = sI.*FI+(iIe),
%  f(iJ) = sJ.*FJ+(iJe),
%  f(iK) = sK.*FK+(iKe),
%  f(iN) = sN.*FN+(iNe),
% where s(I/J/K/N) are sign-vectors. 
% INPUT:
%  G   : grid with field G.cells.indexMap
%  NNC : list og non-neighboring connections in *GLOBAL INDEXING* such that
%        NNC = G.cells.indexMap(NNC_grid). Should be set to [] if there are no 
%        NNCs present, to prevent warning message.
%        
% A later version could concatenate (non-logical) indices, but nnc should
% always come last in case (some) I,J or K -faces are given as nnc 

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

if nargin < 2
    warning('No NNC given ...')
    NNC = [];
end

opt     = struct('neighbors', [], ...
                 'checkConsistency', true);
opt     = merge_options(opt, varargin{:});

nc   = G.cells.num;
dims = G.cartDims;
N    = opt.neighbors;
if isempty(N)
    N = G.faces.neighbors;
end
% switch order s.t. n1 < n2
isPos  = N(:,2) > N(:,1);
N(~isPos,:) = N(~isPos, [2, 1]);
% only consider interior faces
ie = prod(N,2)~=0;

% ----------------- Index to I,J,K - faces ----------------------
% Get global index to interior neighbors (Ng):
%Ng = zeros(G.faces.num, 2);
Ng = zeros(size(N));
Ng(ie,:) = G.cells.indexMap(N(ie,:)); 
ix.iI = and(ie, Ng(:,2)-Ng(:,1) == 1);                 
ix.iJ = and(ie, Ng(:,2)-Ng(:,1) == dims(1));           
ix.iK = and(ie, Ng(:,2)-Ng(:,1) == dims(1)*dims(2));    
% Index to cell-based I,J,K +(plus) faces correspond to first column:
% Need to map back to active grid:
glob2act = zeros(prod(G.cartDims), 1);
glob2act(G.cells.indexMap) = (1:G.cells.num)'; 
ix.iIe = glob2act(Ng(ix.iI,1));
ix.iJe = glob2act(Ng(ix.iJ,1));
ix.iKe = glob2act(Ng(ix.iK,1));
% Also include sign (should always be positive, but just to be sure ...)
ix.sI = 2*isPos(ix.iI)-1;
ix.sJ = 2*isPos(ix.iJ)-1;
ix.sK = 2*isPos(ix.iK)-1;

% ----------------- Index NNC -faces        ----------------------
if ~isempty(NNC)
    % usual mapping from global to active:
    activeIx = zeros(prod(dims),1);
    activeIx(G.cells.indexMap) = (1:G.cells.num)';  
    NNC = activeIx(NNC); % = NNC_grid
    % Switch order s.t. nnc1 < nnc2
    isPosN = NNC(:,2) > NNC(:,1);
    NNC(~isPosN,:) = NNC(~isPosN, [2, 1]);
    nNNC = size(NNC,1);
    % All NNC-cells must correspond to active grid cells:
    assert( all(NNC(:)>0), 'Some NNC-cells are not active in grid' );
    % Produce sparse linear index to NNC-pairs:
    lix = sparse(NNC(:,1) + NNC(:,2)*nc, 1, (1:nNNC)', (nc+1)^2, 1);
    % Find indices to NNCs matching interior neighbors in G:
    [ix.iN, ~, ix.iNe] = find(lix(N(:,1) + N(:,2)*nc));
    % Finally get sign:
    [s1, s2] = deal(2*isPos(ix.iN)-1 , 2*isPosN(ix.iNe)-1);
    ix.sN = s1.*s2;
end
% Perfrom some consistency-checks:
if opt.checkConsistency
    checkConsistency(ix, NNC, ie)
end
end

%% ------------------------------------------------------------------------

function [] = checkConsistency(ix, NNC, ie)
% Check for NNCs not corresponding to grid face in G
if ~isempty(NNC)
    n_fail = size(NNC,1) - numel(ix.iN);
    if n_fail~=0
        warning('Non compatible grid! Found %d NNCs not occuring in grid.', n_fail);
    end
    % Check for non-assigend grid-faces
    ii = ix.iI + ix.iJ + ix.iK;
    ii(ix.iN) = true;
    n_fail = nnz(~ii(ie));
    if n_fail > 0
        warning('Non compatible grid! Found %d interior grid-faces not occuring in output.', n_fail)
    end
end
end



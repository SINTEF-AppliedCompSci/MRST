function Gs = simGridTPFA(G, rock, varargin)
% Construct 'simulation grid' Gs with limited but sufficient info for 
% typical TPFA solvers. Options are 'neighbors', 'porv' and 'depth' 
% referring to grid cells of G. 'neighbors' need not to be compatible 
% with G.faces.   

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

opt = struct('neighbors', [], 'porv', [], 'depth', [], 'actnum', [], 'trans', []);
opt = merge_options(opt, varargin{:});

Gs.type     = 'tpfaGrid';
Gs.cartDims = G.cartDims;
Gs.griddim  = G.griddim;

if isempty(opt.neighbors)
    % faces
    Gs.faces = struct('neighbors', G.faces.neighbors, ...
                      'num',       G.faces.num );  
    % cells
    Gs.cells = struct('facePos',  G.cells.facePos, ...
                      'faces',    G.cells.faces, ...
                      'indexMap', G.cells.indexMap, ...
                      'num',      G.cells.num);
else % use optional neighbor list
    % faces
    N  = opt.neighbors;
    nf = size(N,1);
    Gs.faces = struct('neighbors', N, ...
                      'num',       nf);
    % cells
    if ~isfield(G.cells, 'eMap') || G.cells.eMap == ':'
        nc = G.cells.num;
        indexMap = G.cells.indexMap;
    elseif isnumeric(G.cells.eMap)
        %must have actnum in this case
        assert(~isempty(opt.actnum), 'Numeric ''eMap''-field requires non-empty actnum')
        nc = max(max(G.cells.eMap), max(max(N)));
        indexMap = find(opt.actnum);
    else
        error('Unexpected G.cells.eMap, should be numeric or '':''.');
    end
        
    ix = sortrows([N(:), repmat((1:nf)', [2 1])]);
    ix = ix(ix(:,1)>0,:);
    % Get number of faces per cell
    nfac = zeros(nc, 1);
    [cix, n] = rlencode(ix(:,1));
    nfac(cix) = n;
    facePos = [1; cumsum(nfac)+1];
    faces = ix(:,2);
    Gs.cells = struct('facePos',  facePos, ...
                      'faces',    faces, ...
                      'indexMap', indexMap, ...
                      'num',      nc);
end

% geometry:
if ~isempty(opt.porv)
    Gs.cells.volumes = opt.porv./rock.poro;
    if isfield(rock, 'ntg')
        Gs.cells.volumes = Gs.cells.volumes./rock.ntg;
    end
else
    Gs.cells.volumes = G.cells.volumes;
end

Gs.cells.centroids = nan(G.cells.num, 3);
if isfield(G.cells, 'centroids')
    Gs.cells.centroids(G.cells.eMap,:) = G.cells.centroids;
end
    
if ~isempty(opt.depth)
    Gs.cells.centroids(:,3) = opt.depth;
end

% finally include face areas for some reason ...
Gs.faces.areas = nan(G.faces.num, 1);
if ~isempty(opt.trans)
    Gs.faces.TRANS = opt.trans;
end
end





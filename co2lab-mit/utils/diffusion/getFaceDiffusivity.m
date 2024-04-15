function T = getFaceDiffusivity(G, rock, deck, varargin)
%Compute face diffusivity

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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
    hasDeck = true;
    if nargin == 2 || isempty(deck)
        deck = struct('GRID', struct());
        hasDeck = false;
    end
    
    if mod(numel(varargin), 2) == 1
        varargin = [{deck}, varargin];
        deck = struct('GRID', struct());
        hasDeck = false;
    end
    % Use porosity instead of perm
    local_rock = rock;
    local_rock.perm = rock.poro;
    
    if hasDeck && isfield(deck.GRID, 'COORD')
        [cellCenters,cellFaceCenters]=computeCpGeometry(G,deck.GRID);
        T = computeTrans(G, local_rock,'K_system', 'loc_xyz', ...
                 'cellCenters', cellCenters, ...
                 'cellFaceCenters', cellFaceCenters);
    else
        T = computeTrans(G, local_rock, varargin{:});
    end
    
    % Reduce half-face trans for face transmissibility
    cf = G.cells.faces(:,1);
    nf = G.faces.num;
    T  = 1 ./ accumarray(cf, 1./T, [nf, 1]);
end

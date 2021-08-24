function [G, rock, fluid, deck, state] = setupEGG(varargin)
%Setup the EGG benchmark case for simulation
%
% SYNOPSIS
%   [G, rock, fluid, deck, state] = setupEGG()
%
% PARAMETERS:
%   none
%
% RETURNS:
%   G     - Grid data structure
%   rock  - Petrophysical variables
%   fluid - Fluid structure
%   deck  - Keywords from the input data file
%   state - Initial reservoir state
%
% SEE ALSO:
%   setupSPE1, setupSPE9

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
    opt = struct('useACTNUM', true);
    [opt, extra] = merge_options(opt, varargin{:});
    deck = getDeckEGG(extra{:});
    G = initEclipseGrid(deck);

    G = computeGeometry(G);
    
    rock  = initEclipseRock(deck);
    rock  = compressRock(rock, G.cells.indexMap);

    % Create a special ADI fluid which can produce differentiated fluid
    % properties.
    fluid = initDeckADIFluid(deck);

    % The case includes gravity
    gravity reset on

    % Gravity initialization
    pr   = 400*barsa;
    rz   = G.cells.centroids(1,3);
    dz   = G.cells.centroids(:,3) - rz;
    rhoO    = fluid.bO(400*barsa)*fluid.rhoOS;
    rhoW    = fluid.bW(400*barsa)*fluid.rhoWS;
    rhoMix  = .1*rhoW + .9*rhoO;
    p0   = pr + norm(gravity)*rhoMix*dz;      

    state = initResSol(G, p0, [0.1, .90]);            

end
function out = setupEggFn(n, mode)
% Helper function for the example script ensembleGUIforEgg.m
%
% SYNOPSIS:
%   out = setupEggFn(n, mode)
%
% DESCRIPTION:
%  This function returns the structure required by the ModelEnsemble class
%  setupFn. This contains the necessary structures to run both flow
%  diagnostics and simulations for ensemble member n of the Egg ensemble.
%
% PARAMETERS:
%   n  - the number of the ensemble member required.
%
%   mode - if this is "simulation" (default) then appropriate simulation options
%          will be generated.  
%
% RETURNS:
%     out - structure required by ModelEnsemble().
%    
% EXAMPLE:
%   out = setupEggFn(1);
%
% SEE ALSO:
%   `ensembleGUIForEgg`, `ModelEnsemble`
    
    if nargin < 2
        mode = 'simulation';
    end
    out = struct('model', [], 'state0', [], 'W', [], 'schedule', [], 'simOpts', {{}});

    deck = getDeckEGG('realization', n);
    G = initEclipseGrid(deck);
    G = computeGeometry(G);

    if strcmp(mode, 'simulation')
        [out.state0, out.model, out.schedule, nonlinear] = ...
            initEclipseProblemAD(deck, 'G', G, ...
            'TimestepStrategy', 'none', 'useMex', true);
        out.simOpts  = {'NonLinearSolver', nonlinear};
        out.W = out.schedule.control(1).W;
    else % diagnostics
        rock  = initEclipseRock(deck);
        rock  = compressRock(rock, G.cells.indexMap);
        fluid = initDeckADIFluid(deck);
        out.model = selectModelFromDeck(G, rock, fluid, deck);
        schedule = convertDeckScheduleToMRST(model, deck);
        out.W = schedule.control(1).W;
    end
end

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

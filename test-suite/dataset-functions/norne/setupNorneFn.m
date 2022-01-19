function out = setupNorneFn(n)
% Helper function for setting up a simplified Norne model
%
% SYNOPSIS:
%   out = setupNorneFn(n)
%
% DESCRIPTION:
%  This function returns the structure required by the ModelEnsemble class
%  setupFn. This contains the necessary structures to run both flow
%  diagnostics and simulations for ensemble member n of the Egg ensemble.
%
% PARAMETERS:
%   n  - the number of the ensemble member required.
%
% RETURNS:
%     out - structure required by ModelEnsemble().
%    
% EXAMPLE:
%   out = setupNorneFn(1);
%
% SEE ALSO:
%   `setupNorneRealization`
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

    out = struct('model', [], 'state0', [], 'W', [], 'schedule', [], 'simOpts', {{}});


    % Setup function for Norne ensemble.

    [G, rock, fluid, deck] = setupNorneRealization(n);
    
    deck.RUNSPEC = [];
    deck.PROPS = []; 
    deck.SOLUTION = [];     
    
    model = GenericBlackOilModel(G, rock, fluid, 'water', true, 'oil', true, 'gas', false, 'inputdata', deck);

    W = setupNorneWells(G,rock);    
    yr        = [31 28 31 30 31 30 31 31 30 31 30 31]*day;
    lpyr      = repmat(yr',1,4); lpyr(2,4)=29;
    timesteps = repmat(lpyr(:),3,1);
    schedule  = simpleSchedule(timesteps, 'W', W);
    state0    = initState(G, W, 100*barsa, [0, 1]);
    
    nonlinear = NonLinearSolver();
    
    out.state0 = state0;
    out.model = model;
    out.W = W;
    out.schedule = schedule;

    out.simOpts  = {'NonLinearSolver', nonlinear};
    
    
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
